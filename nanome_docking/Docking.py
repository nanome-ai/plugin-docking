import nanome
from nanome.util import Logs
from ._DockingCalculations import DockingCalculations as Smina
from ._DockingCalculationsAutodock4 import DockingCalculations as Autodock4
from ._DockingCalculationsRhodium import DockingCalculations as Rhodium
from .menus._DockingMenu import DockingMenu
from .menus._DockingMenuRhodium import DockingMenuRhodium
import functools
import sys

from .ComplexUtils import ComplexUtils

__metaclass__ = type
class Docking(nanome.PluginInstance):
    def __init__(self):
        self._menu = None
        self.setting_menu = None
        self._calculations = None
        self._autobox = True

    # Called when Nanome connects to the Plugin, after its instantiation
    def start(self):
        self._menu.build_menu()
        if self._autobox == False:
            self._menu.disable_autobox()
        # Request shallow complex (name, position, orientation), to display them in a list
        self.request_complex_list(self.on_complex_list_received)

    # Called when user clicks on the "Run" button in Nanome
    def on_run(self):
        self.menu.enabled = True
        self.update_menu(self.menu)

    # Called when user click on the "Advanced Settings" button in Nanome
    def on_advanced_settings(self):
        self.setting_menu.enabled = True
        self.setting_menu.index = 1
        self.update_menu(self.setting_menu)

    # Called when a complex is added to the workspace in Nanome
    def on_complex_added(self):
        self.request_complex_list(self.on_complex_list_received)

    # Called when a complex is removed from the workspace in Nanome
    def on_complex_removed(self):
        self.request_complex_list(self.on_complex_list_received)

    def make_plugin_usable(self):
        self._menu.make_plugin_usable()

    def on_complex_list_received(self, complexes):
        self._menu.change_complex_list(complexes)

    def replace_conformer(self, complexes, callback, existing=True):
        for i in range(len(complexes)):
            complex_index = complexes[i].index
            complexes[i] = complexes[i].convert_to_frames()
            complexes[i].index = complex_index

        rerequest_complexes = functools.partial(self.request_complexes, [complex.index for complex in complexes], callback)
        if not existing:
            request_docking_results = functools.partial(self.request_docking_results, callback)
            rerequest_complexes = functools.partial(self.request_complex_list, request_docking_results)

        self.update_structures_deep(complexes, rerequest_complexes)

    def request_docking_results(self, callback, all_complexes):
        docking_results_index = all_complexes[len(all_complexes)-1].index
        self.request_complexes([docking_results_index], callback)

    def set_and_convert_structures(self, has_site, params, complexes):
        # When deep complexes data are received, unpack them and prepare ligand for docking
        receptor = complexes[0]
        self._receptor = receptor
        starting_lig_idx = 1
        site = None

        if has_site:
            site = complexes[1]  #._deep_copy()  # TODO VOODOO
            self._site = site
            ComplexUtils.align_to(site, receptor)
            starting_lig_idx = 2

        # convert ligands to frame representation
        start_docking = lambda ligands: self._calculations.start_docking(receptor, ligands, site, **params)
        self.replace_conformer(complexes[starting_lig_idx:], start_docking)

    def run_docking(self, receptor, ligands, site, params):
        # Change the plugin to be "unusable"
        if self._menu._run_button.unusable == True:
            return
        self._menu.make_plugin_usable(False)
        #self._menu.show_loading(True)

        # Request complexes to Nanome in this order: [receptor, site (if any), ligand, ligand,...]
        request_list = [receptor.index]
        if site != None:
            request_list.append(site.index)
        request_list += [x.index for x in ligands]

        setup_structures = functools.partial(self.set_and_convert_structures, site != None, params)
        self.request_complexes(request_list, setup_structures)
        #self._menu.show_loading(False)

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


def main():
    name = None
    cl = None
    for arg in sys.argv:
        if arg == "smina":
            name = "Smina"
            cl = SminaDocking
        elif arg == "autodock4":
            name = "Autodock 4"
            cl = Autodock4Docking
        elif arg == "rhodium":
            name = "Rhodium"
            cl = RhodiumDocking
    if name == None:
        Logs.error("Please pass the docking software to use as an argument: smina|autodock4|rhodium")
        sys.exit(1)

    # Create the plugin, register Docking as the class to instantiate, and start listening
    plugin = nanome.Plugin(name + " Docking", "Run docking using " + name + ". Lets user choose the receptor, ligands, and diverse options", "Docking", True)
    plugin.set_plugin_class(cl)
    plugin.run()

if __name__ == "__main__":
    main()
