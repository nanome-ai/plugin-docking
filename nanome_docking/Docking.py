import nanome
from nanome.util import async_callback, ComplexUtils

from nanome_docking.smina.calculations import DockingCalculations as Smina
from nanome_docking.autodock4.calculations import DockingCalculations as Autodock4
from nanome_docking.rhodium.calculations import DockingCalculations as Rhodium
from nanome_docking.menus.DockingMenu import DockingMenu, SettingsMenu
from nanome_docking.menus.DockingMenuRhodium import DockingMenuRhodium


__metaclass__ = type


class Docking(nanome.AsyncPluginInstance):

    def __init__(self):
        self.menu = None
        self._calculations = None

    @async_callback
    async def start(self):
        self.settings_menu = SettingsMenu(self)
        self.menu.build_menu()
        # Request shallow complex (name, position, orientation), to display them in a list
        complexes = await self.request_complex_list()
        self.menu.change_complex_list(complexes)

    def on_run(self):
        # Called when user clicks on the "Run" button in Nanome
        self.menu.enable()

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
        self._calculations.start_docking(receptor, ligands, site, **params)

    def add_result_to_workspace(self, results, align=False):
        for complex in results:
            complex.position = self._receptor.position
            complex.rotation = self._receptor.rotation

            if align and hasattr(self, '_site'):
                ComplexUtils.align_to(complex, self._site)
                complex.boxed = True

        self.update_structures_deep(results)


class SminaDocking(Docking):

    def __init__(self):
        super(SminaDocking, self).__init__()
        self.menu = DockingMenu(self)
        self._calculations = Smina(self)


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
