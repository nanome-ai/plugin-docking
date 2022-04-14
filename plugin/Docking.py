import nanome
from nanome.util import async_callback, ComplexUtils
import os
import re
import tempfile

from nanome.util.enums import NotificationTypes
from nanome.util import Logs

from plugin.smina.calculations import DockingCalculations as Smina
from plugin.autodock4.calculations import DockingCalculations as Autodock4
from plugin.menus.DockingMenu import DockingMenu, SettingsMenu

__metaclass__ = type


PDBOPTIONS = nanome.api.structure.Complex.io.PDBSaveOptions()
PDBOPTIONS.write_bonds = True

# Default Timeout should be 5 minutes per frame
DEFAULT_TIMEOUT = 300
TIMEOUT_PER_FRAME = int(os.environ.get('TIMEOUT_PER_FRAME', DEFAULT_TIMEOUT))


class Docking(nanome.AsyncPluginInstance):

    def __init__(self):
        super().__init__()
        self.menu = DockingMenu(self)
        self.settings_menu = SettingsMenu(self)
        self.docked_complex_indices = []

    def start(self):
        self.menu.build_menu()

    @async_callback
    async def on_run(self):
        # Called when user clicks on the "Run" button in Nanome
        self.menu.enable()
        complexes = await self.request_complex_list()
        self.menu.change_complex_list(complexes)

    def on_advanced_settings(self):
        # Called when user click on the "Advanced Settings" button in Nanome
        self.settings_menu.enable()

    @async_callback
    async def on_complex_added(self):
        # Called when a complex is added to the workspace in Nanome
        complexes = await self.request_complex_list()
        self.menu.change_complex_list(complexes)

    @async_callback
    async def on_complex_removed(self):
        # Called when a complex is removed from the workspace in Nanome
        complexes = await self.request_complex_list()
        self.menu.change_complex_list(complexes)
        # If docked complex has been deleted, remove index from list
        comp_indices = [cmp.index for cmp in complexes]
        self.docked_complex_indices = [x for x in self.docked_complex_indices if x in comp_indices]

    async def run_docking(self, receptor, ligands, site, params):
        # Request complexes to Nanome in this order: [receptor, <site>, ligand, ligand,...]
        # site not always required.
        complex_indices = [receptor.index]
        if site:
            complex_indices += [site.index]
        complex_indices += [x.index for x in ligands]
        complexes = await self.request_complexes(complex_indices)
        receptor = complexes[0]

        if site:
            site = complexes[1]
            ligands = complexes[2:]
        else:
            ligands = complexes[1:]

        ComplexUtils.convert_to_frames(ligands)

        # Get advanced_settings.
        advanced_settings = self.settings_menu.get_settings()
        params.update(advanced_settings)

        with tempfile.TemporaryDirectory() as temp_dir:
            # Convert input complexes into PDBs.
            receptor_pdb = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb", dir=temp_dir)
            site_pdb = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb", dir=temp_dir)
            receptor.io.to_pdb(receptor_pdb.name, PDBOPTIONS)
            site.io.to_pdb(site_pdb.name, PDBOPTIONS)

            ligand_pdbs = []
            for lig in ligands:
                cleaned_name = lig.full_name.replace(' ', '_')
                ligand_pdb = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb", dir=temp_dir, prefix=f'{cleaned_name}_')
                ComplexUtils.align_to(lig, receptor)
                lig.io.to_pdb(ligand_pdb.name, PDBOPTIONS)
                ligand_pdbs.append(ligand_pdb)

            self.send_notification(NotificationTypes.message, "Docking started")
            output_complexes = []
            frame_count = 0
            for lig in ligands:
                frame_count += sum(1 for _ in lig.molecules)
            timeout = TIMEOUT_PER_FRAME * frame_count
            Logs.message(f'Running Docking on {len(ligands)} ligands, containing a total of {frame_count} frames'.format(frame_count, timeout))
            try:
                output_sdfs = await self._calculations.start_docking(
                    receptor_pdb, ligand_pdbs, site_pdb, temp_dir, timeout=timeout, **params)
            except TimeoutError:
                message = "Docking run timed out"
                self.send_notification(NotificationTypes.error, message)
                Logs.error(message)
                return

            for ligand, result in zip(ligands, output_sdfs):
                docked_complex = nanome.structure.Complex.io.from_sdf(path=result.name)
                if len(list(docked_complex.molecules)) == 0:
                    msg = "Docking run returned 0 results."
                    Logs.warning(msg)
                    self.send_notification(NotificationTypes.warning, msg)
                    return

                docked_complex.full_name = f'{ligand.full_name} (Docked)'
                docked_complex = docked_complex.convert_to_frames()
                # fix metadata sorting
                if hasattr(self, 'set_scores'):
                    for molecule in docked_complex.molecules:
                        self.set_scores(molecule)

                show_atom_labels = params.get('visual_scores', False)
                if hasattr(self, 'visualize_scores'):
                    self.visualize_scores(docked_complex, show_atom_labels=show_atom_labels)

                docked_complex.set_current_frame(0)
                docked_complex.visible = True
                docked_complex.locked = True
                output_complexes.append(docked_complex)

        # hide ligands
        for ligand in ligands:
            ligand.visible = False
            ComplexUtils.reset_transform(ligand)
        self.update_structures_shallow(ligands)

        # Add docked complexes to workspace.
        await self.add_result_to_workspace(output_complexes, receptor, site)
        self.send_notification(NotificationTypes.success, "Docking finished")
        return output_complexes

    async def add_result_to_workspace(self, results, receptor, site):
        for comp in results:
            comp.position = receptor.position
            comp.rotation = receptor.rotation
            ComplexUtils.align_to(comp, site)
            comp.boxed = True
        created_complexes = await self.add_to_workspace(results)
        # add_to_workspace doesn't set correct current frame, so set it back to 0.
        # Also need to manually reset position and rotation
        for c1, c2 in zip(created_complexes, results):
            c1.position = c2.position
            c1.rotation = c2.rotation
            c1.set_current_frame(0)
        self.update_structures_shallow(created_complexes)
        indices = [cmp.index for cmp in created_complexes]
        self.docked_complex_indices.extend(indices)

    def enable_loading_bar(self, enabled=True):
        self.menu.enable_loading_bar(enabled)

    def update_loading_bar(self, current, total):
        self.menu.update_loading_bar(current, total)

    def update_run_btn_text(self, new_text):
        self.menu.update_run_btn_text(new_text)

    async def toggle_atom_labels(self, enabled: bool):
        docked_complexes = await self.request_complexes(self.docked_complex_indices)
        for comp in docked_complexes:
            for atom in comp.atoms:
                atom.labeled = bool(enabled and atom.label_text)
        await self.update_structures_deep(docked_complexes)


class SminaDocking(Docking):

    def __init__(self):
        super(SminaDocking, self).__init__()
        self.menu = DockingMenu(self)
        self._calculations = Smina(self)

    def set_scores(self, molecule):
        """Clean and Parse score information for provided molecule."""
        molecule.min_atom_score = float('inf')
        molecule.max_atom_score = float('-inf')

        num_rgx = r'(-?[\d.]+(?:e[+-]\d+)?)'
        pattern = re.compile('<{},{},{}> {} {} {} {} {}'.format(*([num_rgx] * 8)), re.U)
        for associated in molecule.associateds:
            # make the labels pretty :)
            associated['Minimized Affinity'] = associated.pop('> <minimizedAffinity>')
            associated['Atomic Interaction Terms'] = associated.pop('> <atomic_interaction_terms>')

            pose_score = associated['Minimized Affinity']
            for residue in molecule.residues:
                residue.label_text = pose_score + " kcal/mol"
                residue.labeled = True
            interaction_terms = associated['Atomic Interaction Terms']
            interaction_values = re.findall(pattern, interaction_terms)
            for i, atom in enumerate(molecule.atoms):
                if i < len(interaction_values) - 1:
                    # Logs.debug("interaction values for atom " + str(i) + ": " + str(interaction_values[i]))
                    atom.score = float(interaction_values[i][5])
                    molecule.min_atom_score = min(atom.score, molecule.min_atom_score)
                    molecule.max_atom_score = max(atom.score, molecule.max_atom_score)

    def visualize_scores(self, ligand_complex, show_atom_labels=False):
        for molecule in ligand_complex.molecules:
            for atom in molecule.atoms:
                if hasattr(atom, "score") and atom.score != 0.0:
                    atom.label_text = self._truncate(atom.score, 3)
                    atom.labeled = show_atom_labels

    def _truncate(self, f, n):
        """Truncates/pads a float f to n decimal places without rounding."""
        s = '{}'.format(f)
        if 'e' in s or 'E' in s:
            return '{0:.{1}f}'.format(f, n)
        i, p, d = s.partition('.')
        return '.'.join([i, (d + '0' * n)[:n]])


class Autodock4Docking(Docking):

    def __init__(self):
        super(Autodock4Docking, self).__init__()
        self.menu = DockingMenu(self)
        self._calculations = Autodock4(self)

    def set_scores(self, molecule):
        for associated in molecule.associateds:
            associated.pop('>  <MODEL>')
            associated.pop('>  <TORSDO>')
            remark = associated.pop('>  <REMARK>')

            split_remark = remark.split()
            associated['CONF_DEPENDENT'] = split_remark[3]
            associated['TOTAL_SCORE'] = split_remark[4]
            associated['INTER + INTRA'] = split_remark[8]
            associated['INTER'] = split_remark[10]
            associated['INTRA'] = split_remark[12]
            associated['CONF_INDEPENDENT'] = split_remark[14]
            associated['UNBOUND'] = split_remark[16]
