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

    def initialize(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self._receptor_input = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb", dir=self.temp_dir.name)
        self._ligands_input = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb", dir=self.temp_dir.name)
        self._site_input = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb", dir=self.temp_dir.name)
        self._docking_output = tempfile.NamedTemporaryFile(delete=False, prefix="output", suffix=".sdf", dir=self.temp_dir.name)
        self._ligand_output = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb", dir=self.temp_dir.name)
        self._log_file = tempfile.NamedTemporaryFile(delete=False, dir=self.temp_dir.name)

    async def start_docking(self, receptor, ligands, site, exhaustiveness, modes, align, replace, scoring, visual_scores, autobox):
        self._exhaustiveness = exhaustiveness
        self._modes = modes
        self._receptor = receptor
        self._ligands = ligands
        self._combined_ligands = ComplexUtils.combine_ligands(receptor, ligands)
        self._site = site
        self._align = align
        self._replace = replace
        self._scoring = scoring
        self._visual_scores = visual_scores
        self._autobox = autobox

        # Start docking process
        self._write_structures_to_file()
        await self._start_docking()

    def _write_structures_to_file(self):
        self.initialize()
        # Save all input files
        self._receptor.io.to_pdb(self._receptor_input.name, PDBOPTIONS)
        Logs.debug("Saved PDB", self._receptor_input.name)
        self._combined_ligands.io.to_pdb(self._ligands_input.name, PDBOPTIONS)
        Logs.debug("Saved PDB", self._ligands_input.name)
        self._site.io.to_pdb(self._site_input.name, PDBOPTIONS)
        Logs.debug("Saved PDB", self._site_input.name)

    async def _start_docking(self):
        smina_args = [
            '-r', self._receptor_input.name,
            '-l', self._ligands_input.name,
            '--autobox_ligand', self._site_input.name]

        if self._scoring:
            smina_args += ['--score_only', '--out', self._ligand_output.name]
        else:
            smina_args += [
                '--out', self._docking_output.name,
                '--log', self._log_file.name,
                '--exhaustiveness', str(self._exhaustiveness),
                '--num_modes', str(self._modes),
                '--autobox_add', str(self._autobox),
                '--atom_term_data']

        Logs.debug("Run SMINA")
        self._start_timer = timer()

        p = Process(SMINA_PATH, smina_args)
        p.on_done = self._docking_finished
        p.on_output = self.update_loading_bar
        self.plugin.send_notification(NotificationTypes.message, "Docking started")
        await p.start()

    def update_loading_bar(self, output):
        # Smina program outputs a loading bar to stdout using asterisks.
        # Update UI loading bar based on number of asterisks outputted
        total_stars = 51
        output_str = output.decode()
        star_count = output_str.count('*')
        if star_count:
            self.loading_bar_counter += star_count
            self.plugin.update_loading_bar(self.loading_bar_counter, total_stars)

    def _docking_finished(self, return_code):
        if return_code != 0:
            message = 'Docking error, check plugin'
            self.plugin.send_notification(NotificationTypes.error, message)
            Logs.error(message)
            return

        end = timer()
        Logs.debug("Docking Finished in", end - self._start_timer, "seconds")

        # hide ligands
        for ligand in self._ligands:
            ligand.visible = False
            ComplexUtils.reset_transform(ligand)
        self.plugin.update_structures_shallow(self._ligands)

        docking_results = nanome.structure.Complex.io.from_sdf(path=self._docking_output.name)
        ComplexUtils.convert_to_frames([docking_results])

        # fix metadata sorting
        docking_results._remarks['Minimized Affinity'] = ''

        Logs.debug("Read SDF", self._docking_output.name)
        for molecule in docking_results.molecules:
            self._set_scores(molecule)

        docking_results.set_current_frame(0)

        if self._visual_scores:
            self._visualize_scores(docking_results)

        if not self._scoring:
            if len(self._combined_ligands.names) > 1:
                docking_results.name += "Docking Results"
            elif len(self._combined_ligands.names) == 1:
                docking_results.name = self._combined_ligands.names[0] + " (Docked)"
            docking_results.visible = True
            docking_results.locked = True

        if self._scoring:
            Logs.debug("Display scoring result")
            self.plugin.display_scoring_result(docking_results)
        else:
            Logs.debug("Update workspace")
            Logs.debug(f'** result index {docking_results.index}')
            ComplexUtils.convert_to_conformers([docking_results])
            self.plugin.add_result_to_workspace([docking_results], self._align)

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
