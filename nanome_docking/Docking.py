import argparse

import nanome
from nanome.util import async_callback, ComplexUtils

from nanome_docking.smina.calculations import DockingCalculations as Smina
from nanome_docking.autodock4.calculations import DockingCalculations as Autodock4
from nanome_docking.rhodium.calculations import DockingCalculations as Rhodium
from nanome_docking.menus.DockingMenu import DockingMenu
from nanome_docking.menus.DockingMenuRhodium import DockingMenuRhodium


__metaclass__ = type


class Docking(nanome.AsyncPluginInstance):
    def __init__(self):
        self._menu = None
        self.setting_menu = None
        self._calculations = None
        self._autobox = True

    @async_callback
    async def start(self):
        # Called when Nanome connects to the Plugin, after its instantiation
        self._menu.build_menu()
        if self._autobox is False:
            self._menu.disable_autobox()
        # Request shallow complex (name, position, orientation), to display them in a list
        complexes = await self.request_complex_list()
        self._menu.change_complex_list(complexes)

    def on_run(self):
        # Called when user clicks on the "Run" button in Nanome
        self.menu.enabled = True
        self.update_menu(self.menu)

    def on_advanced_settings(self):
        # Called when user click on the "Advanced Settings" button in Nanome
        self.setting_menu.enabled = True
        self.setting_menu.index = 1
        self.update_menu(self.setting_menu)

    @async_callback
    async def on_complex_added(self):
        # Called when a complex is added to the workspace in Nanome
        complexes = await self.request_complex_list()
        self._menu.change_complex_list(complexes)

    @async_callback
    async def on_complex_removed(self):
        # Called when a complex is removed from the workspace in Nanome
        complexes = await self.request_complex_list()
        self._menu.change_complex_list(complexes)

    def make_plugin_usable(self):
        self._menu.make_plugin_usable()

    def set_and_convert_structures(self, has_site, params, complexes):
        # When deep complexes data are received, unpack them and prepare ligand for docking
        receptor = complexes[0]
        self._receptor = receptor
        starting_lig_idx = 1
        site = None

        if has_site:
            site = complexes[1]
            self._site = site
            ComplexUtils.align_to(site, receptor)
            starting_lig_idx = 2

        # convert ligands to frame representation
        ligands = complexes[starting_lig_idx:]
        ComplexUtils.convert_to_frames(ligands)
        self._calculations.start_docking(receptor, ligands, site, **params)

    async def run_docking(self, receptor, ligands, site, params):
        # Change the plugin to be "unusable"
        if self._menu._run_button.unusable is True:
            return
        self._menu.make_plugin_usable(False)
        # self._menu.show_loading(True)

        # Request complexes to Nanome in this order: [receptor, site (if any), ligand, ligand,...]
        request_list = [receptor.index]
        if site is not None:
            request_list.append(site.index)
        request_list += [x.index for x in ligands]

        complexes = await self.request_complexes(request_list)
        has_site = site is not None
        self.set_and_convert_structures(has_site, params, complexes)
        # self._menu.show_loading(False)

    # Called every update tick of the Plugin
    def update(self):
        self._calculations.update()

    def add_result_to_workspace(self, results, align=False):
        for complex in results:
            complex.position = self._receptor.position
            complex.rotation = self._receptor.rotation

            if align:
                ComplexUtils.align_to(complex, self._site)
                complex.boxed = True

        self.update_structures_deep(results)

    def display_scoring_result(self, result):
        self._menu.display_scoring_result(result)


class SminaDocking(Docking):
    def __init__(self):
        super(SminaDocking, self).__init__()
        self._calculations = Smina(self)
        self._menu = DockingMenu(self)


class Autodock4Docking(Docking):
    def __init__(self):
        super(Autodock4Docking, self).__init__()
        self._calculations = Autodock4(self)
        self._menu = DockingMenu(self)
        self._autobox = False


class RhodiumDocking(Docking):
    def __init__(self):
        super(RhodiumDocking, self).__init__()
        self._calculations = Rhodium(self)
        self._menu = DockingMenuRhodium(self)
