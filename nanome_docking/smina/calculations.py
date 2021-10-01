import re
import sys
import nanome
import os
import tempfile
import subprocess
from timeit import default_timer as timer

from nanome.util import ComplexUtils, Logs
from nanome.util.enums import NotificationTypes

PDBOPTIONS = nanome.api.structure.Complex.io.PDBSaveOptions()
PDBOPTIONS.write_bonds = True

SMINA_PATH = os.path.join(os.getcwd(), 'nanome_docking', 'smina', 'smina_binary')


class DockingCalculations():

    def __init__(self, plugin):
        self.plugin = plugin
        self.requires_site = True
        self.loading_bar_counter = 0

    async def start_docking(self, receptor, ligands, site, exhaustiveness, modes, align, replace, scoring, visual_scores, autobox):
        # Start docking process
        with tempfile.TemporaryDirectory() as temp_dir:
            receptor_input = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb", dir=temp_dir)
            site_input = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb", dir=temp_dir)
            docking_output = tempfile.NamedTemporaryFile(delete=False, prefix="output", suffix=".sdf", dir=temp_dir)
            log_file = tempfile.NamedTemporaryFile(delete=False, dir=temp_dir)

            receptor.io.to_pdb(receptor_input.name, PDBOPTIONS)
            site.io.to_pdb(site_input.name, PDBOPTIONS)

            # Save all input files
            docking_outputs = []
            _start_timer = timer()
            self.plugin.send_notification(NotificationTypes.message, "Docking started")
            for lig in ligands:
                ligands_input = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb", dir=temp_dir)
                ComplexUtils.align_to(lig, receptor)
                lig.io.to_pdb(ligands_input.name, PDBOPTIONS)

                ligand_output = tempfile.NamedTemporaryFile(delete=False, suffix=".sdf", dir=temp_dir)
                smina_args = [
                    '-r', receptor_input.name,
                    '-l', ligands_input.name,
                    '--autobox_ligand', site_input.name,
                    '--out', ligand_output.name,
                    '--log', log_file.name,
                    '--exhaustiveness', str(exhaustiveness),
                    '--num_modes', str(modes),
                    '--autobox_add', str(autobox),
                    '--atom_term_data'
                ]

                Logs.debug("Run SMINA")

                cmd = [SMINA_PATH, *smina_args]
                process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
                self.handle_loading_bar(process, len(ligands))
                docking_outputs.append(ligand_output)
            end = timer()
            Logs.debug("Docking Finished in", end - _start_timer, "seconds")
            
            # hide ligands
            for ligand in ligands:
                ligand.visible = False
                ComplexUtils.reset_transform(ligand)
            self.plugin.update_structures_shallow(ligands)

            output_complexes = []
            for ligand, result in zip(ligands, docking_outputs):
                docked_complex = nanome.structure.Complex.io.from_sdf(path=result.name)
                docked_complex.full_name = f'{ligand.full_name} (Docked)'
                ComplexUtils.convert_to_frames([docked_complex])
                # fix metadata sorting
                docked_complex._remarks['Minimized Affinity'] = ''

                Logs.debug("Read SDF", docking_output.name)
                for molecule in docked_complex.molecules:
                    self._set_scores(molecule)

                docked_complex.set_current_frame(0)
                docked_complex.visible = True
                docked_complex.locked = True
                output_complexes.append(docked_complex)

            ComplexUtils.convert_to_conformers(output_complexes)
            self.plugin.add_result_to_workspace(output_complexes, align)
            self.plugin.send_notification(NotificationTypes.success, "Docking finished")

    def handle_loading_bar(self, process, ligand_count):
        """Render loading bar from stdout on the menu.

        stdout has a loading bar of asterisks. Every asterisk represents about 2% completed
        """
        stars_per_complex = 51
        total_stars = stars_per_complex * ligand_count

        for c in iter(lambda: process.stdout.read(1), b''):
            if c.decode() == '*':
                self.loading_bar_counter += 1
                self.plugin.update_loading_bar(self.loading_bar_counter, total_stars)
            sys.stdout.buffer.write(c)

    def _set_scores(self, molecule):
        molecule.min_atom_score = float('inf')
        molecule.max_atom_score = float('-inf')

        num_rgx = '(-?[\d.]+(?:e[+-]\d+)?)'
        pattern = re.compile('<{},{},{}> {} {} {} {} {}'.format(*([num_rgx] * 8)), re.U)
        for associated in molecule.associateds:
            # make the labels pretty :)
            associated['Minimized Affinity'] = associated['> <minimizedAffinity>']
            associated['Atomic Interaction Terms'] = associated['> <atomic_interaction_terms>']
            del associated['> <minimizedAffinity>']
            del associated['> <atomic_interaction_terms>']

            pose_score = associated['Minimized Affinity']
            for residue in molecule.residues:
                residue.label_text = pose_score + " kcal/mol"
                residue.labeled = True
            interaction_terms = associated['Atomic Interaction Terms']
            interaction_values = re.findall(pattern, interaction_terms)
            for i, atom in enumerate(molecule.atoms):
                if i < len(interaction_values) - 1:
                    Logs.debug("interaction values for atom " + str(i) + ": " + str(interaction_values[i]))
                    atom.score = float(interaction_values[i][5])
                    molecule.min_atom_score = min(atom.score, molecule.min_atom_score)
                    molecule.max_atom_score = max(atom.score, molecule.max_atom_score)

    def _visualize_scores(self, ligand_complex):
        for molecule in ligand_complex.molecules:
            for atom in molecule.atoms:
                if hasattr(atom, "score"):
                    atom.atom_mode = nanome.api.structure.Atom.AtomRenderingMode.Point
                    denominator = -molecule.min_atom_score if atom.score < 0 else molecule.max_atom_score
                    norm_score = atom.score / denominator
                    atom.atom_scale = norm_score * 1.5 + 0.1
                    atom.label_text = self._truncate(atom.score, 3)
                    # atom.labeled = True

    def _truncate(self, f, n):
        '''Truncates/pads a float f to n decimal places without rounding'''
        s = '{}'.format(f)
        if 'e' in s or 'E' in s:
            return '{0:.{1}f}'.format(f, n)
        i, p, d = s.partition('.')
        return '.'.join([i, (d + '0' * n)[:n]])
