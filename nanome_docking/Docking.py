import nanome
from nanome.util import async_callback, ComplexUtils
import re
import tempfile

from timeit import default_timer as timer
from nanome.util.enums import NotificationTypes
from nanome.util import Logs

from nanome_docking.smina.calculations import DockingCalculations as Smina
from nanome_docking.autodock4.calculations import DockingCalculations as Autodock4
from nanome_docking.rhodium.calculations import DockingCalculations as Rhodium
from nanome_docking.menus.DockingMenu import DockingMenu, SettingsMenu
from nanome_docking.menus.DockingMenuRhodium import DockingMenuRhodium

__metaclass__ = type


class Docking(nanome.AsyncPluginInstance):

    def __init__(self):
        self.menu = DockingMenu(self)
        self.settings_menu = SettingsMenu(self)

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

    async def run_docking(self, receptor, ligands, site, params):
        # Request complexes to Nanome in this order: [receptor, <site>, ligand, ligand,...]
        # site not always required.
        complex_indices = [receptor.index]
        if site:
            complex_indices += [site.index]
        complex_indices += [x.index for x in ligands]

        complexes = await self.request_complexes(complex_indices)
        receptor = complexes[0]
        self._receptor = receptor

        if site:
            site = complexes[1]
            ligands = complexes[2:]
        else:
            ligands = complexes[1:]

        ComplexUtils.convert_to_frames(ligands)

        start_timer = timer()
        self.send_notification(NotificationTypes.message, "Docking started")
        output_complexes = []
        with tempfile.TemporaryDirectory() as temp_dir:
            output_sdfs = await self._calculations.start_docking(receptor, ligands, site, temp_dir, **params)
            end = timer()
            Logs.debug("Docking Finished in", end - start_timer, "seconds")

            for ligand, result in zip(ligands, output_sdfs):
                docked_complex = nanome.structure.Complex.io.from_sdf(path=result.name)
                docked_complex.full_name = f'{ligand.full_name} (Docked)'
                ComplexUtils.convert_to_frames([docked_complex])
                # fix metadata sorting
                docked_complex._remarks['Minimized Affinity'] = ''
                if hasattr(self, '_set_scores'):
                    for molecule in docked_complex.molecules:
                        self._set_scores(molecule)

                docked_complex.set_current_frame(0)
                docked_complex.visible = True
                docked_complex.locked = True
                output_complexes.append(docked_complex)

        # hide ligands
        for ligand in ligands:
            ligand.visible = False
            ComplexUtils.reset_transform(ligand)
        self.update_structures_shallow(ligands)
    
        ComplexUtils.convert_to_conformers(output_complexes)
        self.add_result_to_workspace(output_complexes, receptor, site)
        self.send_notification(NotificationTypes.success, "Docking finished")

    def add_result_to_workspace(self, results, receptor, site):
        for complex in results:
            complex.position = receptor.position
            complex.rotation = receptor.rotation
            ComplexUtils.align_to(complex, site)
            complex.boxed = True
        self.update_structures_deep(results)

    def enable_loading_bar(self, enabled=True):
        self.menu.enable_loading_bar(enabled)

    def update_loading_bar(self, current, total):
        self.menu.update_loading_bar(current, total)


class SminaDocking(Docking):

    def __init__(self):
        super(SminaDocking, self).__init__()
        self.menu = DockingMenu(self)
        self._calculations = Smina(self)

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

class Autodock4Docking(Docking):

    def __init__(self):
        super(Autodock4Docking, self).__init__()
        self.menu = DockingMenu(self)
        self._calculations = Autodock4(self)


class RhodiumDocking(Docking):

    def __init__(self):
        super(RhodiumDocking, self).__init__()
        self.menu = DockingMenuRhodium(self)
        self._calculations = Rhodium(self)
