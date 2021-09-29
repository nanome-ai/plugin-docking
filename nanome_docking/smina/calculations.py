import re
import nanome
import os
import tempfile
from timeit import default_timer as timer

from nanome.util import ComplexUtils, Logs, Process
from nanome.util.asyncio import async_callback
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
        _exhaustiveness = exhaustiveness
        _modes = modes
        _receptor = receptor
        _ligands = ligands
        _combined_ligands = ComplexUtils.combine_ligands(receptor, ligands)
        _site = site
        _align = align
        _replace = replace
        _scoring = scoring
        _visual_scores = visual_scores
        _autobox = autobox

        # Start docking process
        temp_dir = tempfile.TemporaryDirectory()
        _receptor_input = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb", dir=temp_dir.name)
        _ligands_input = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb", dir=temp_dir.name)
        _site_input = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb", dir=temp_dir.name)
        _docking_output = tempfile.NamedTemporaryFile(delete=False, prefix="output", suffix=".sdf", dir=temp_dir.name)
        _ligand_output = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb", dir=temp_dir.name)
        _log_file = tempfile.NamedTemporaryFile(delete=False, dir=temp_dir.name)

        # Save all input files
        _receptor.io.to_pdb(_receptor_input.name, PDBOPTIONS)
        Logs.debug("Saved PDB", _receptor_input.name)
        _combined_ligands.io.to_pdb(_ligands_input.name, PDBOPTIONS)
        Logs.debug("Saved PDB", _ligands_input.name)
        _site.io.to_pdb(_site_input.name, PDBOPTIONS)
        Logs.debug("Saved PDB", _site_input.name)
        
        smina_args = [
            '-r', _receptor_input.name,
            '-l', _ligands_input.name,
            '--autobox_ligand', _site_input.name,
            '--out', _docking_output.name,
            '--log', _log_file.name,
            '--exhaustiveness', str(_exhaustiveness),
            '--num_modes', str(_modes),
            '--autobox_add', str(_autobox),
            '--atom_term_data'
        ]

        Logs.debug("Run SMINA")
        _start_timer = timer()

        import subprocess
        cmd = [SMINA_PATH, *smina_args]
        self.plugin.send_notification(NotificationTypes.message, "Docking started")
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        # stdout has a loading bar of asterisks. Every asterisk represents about 2% completed
        # update loading bar on menu accordingly
        star_count = 0
        total_stars = 51

        import sys
        for c in iter(lambda: process.stdout.read(1), b''):
            if c.decode() == '*':
                star_count += 1
                self.plugin.update_loading_bar(star_count, total_stars)
            sys.stdout.buffer.write(c)

        end = timer()
        Logs.debug("Docking Finished in", end - _start_timer, "seconds")

        # hide ligands
        for ligand in _ligands:
            ligand.visible = False
            ComplexUtils.reset_transform(ligand)
        self.plugin.update_structures_shallow(_ligands)

        docking_results = nanome.structure.Complex.io.from_sdf(path=_docking_output.name)
        ComplexUtils.convert_to_frames([docking_results])

        # fix metadata sorting
        docking_results._remarks['Minimized Affinity'] = ''

        Logs.debug("Read SDF", _docking_output.name)
        for molecule in docking_results.molecules:
            self._set_scores(molecule)

        docking_results.set_current_frame(0)

        if len(_combined_ligands.names) > 1:
            docking_results.name += "Docking Results"
        elif len(_combined_ligands.names) == 1:
            docking_results.name = _combined_ligands.names[0] + " (Docked)"
        docking_results.visible = True
        docking_results.locked = True

        Logs.debug("Update workspace")
        Logs.debug(f'** result index {docking_results.index}')
        ComplexUtils.convert_to_conformers([docking_results])
        self.plugin.add_result_to_workspace([docking_results], _align)
        self.plugin.send_notification(NotificationTypes.success, "Docking finished")

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
            atom_count = len(list(molecule.atoms))
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
