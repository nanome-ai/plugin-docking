import nanome
from _DockingCalculations import DockingCalculations as Smina
from _DockingCalculationsAutodock4 import DockingCalculations as Autodock4
from _DockingMenu import DockingMenu
import sys

__metaclass__ = type
class Docking(nanome.PluginInstance):
    def __init__(self):
        self._menu = DockingMenu(self)
        self._calculations = None

    # Function called when Nanome connects to the Plugin, after its instantiation
    def start(self):
        self._menu.build_menu()

    # Function called when user clicks on the "Run" button in Nanome
    def on_run(self):
        menu = self._menu
        if menu._selected_receptor == None or menu._selected_site == None or menu._selected_ligands == []:
            self.open_menu()
        else:
            self.run_docking(menu._selected_receptor[1], menu._selected_ligands, menu._selected_site[1])

    def on_advanced_settings(self):
        nanome.util.Logs.debug("Advanced Settings")
        self.open_menu()

    def request_refresh(self):
        self.request_complex_list(self.on_complex_list_received)
        nanome.util.Logs.debug("Complex list requested")

    def open_menu(self):
        menu = nanome.ui.Menu.get_plugin_menu()
        menu.enabled = True
        self.update_menu(menu)

    def make_plugin_usable(self):
        self._menu.make_plugin_usable()

    # Function called when Nanome returns the complex list after a request
    def on_complex_list_received(self, complexes):
        nanome.util.Logs.debug("Complex list received:", complexes)
        self._menu.change_complex_list(complexes)

    def run_docking(self, receptor, ligand_list, site):
        if self._menu._run_button.unusable == True:
            return
        self._menu.make_plugin_usable(False)
        request_list = [receptor.index, site.index] + [x.index for x in ligand_list]
        self.request_complexes(request_list, self.on_complexes_received)

    def on_complexes_received(self, complexes):
        receptor = complexes[0]
        site = complexes[1]
        ligands = nanome.structure.Complex()
        for ligand in complexes[2:]:
            for molecule in ligand.molecules:
                ligands.add_molecule(molecule)

        self._calculations.start_docking(receptor, ligands, site, self._menu._exhaustiveness, self._menu._modes, self._menu._align, self._menu._replace, self._menu._scoring_only, self._menu._autobox_size)

    # Function called every update tick of the Plugin
    def update(self):
        self._calculations.update()

    def add_result_to_workspace(self, result):
        self.add_to_workspace(result)

    def display_scoring_result(self, result):
        self._menu.display_scoring_result(result)

class SminaDocking(Docking):
    def __init__(self):
        super(SminaDocking, self).__init__()
        self._calculations = Smina(self)

class Autodock4Docking(Docking):
    def __init__(self):
        super(Autodock4Docking, self).__init__()
        self._calculations = Autodock4(self)

if __name__ == "__main__":
    name = None
    cl = None
    for arg in sys.argv:
        if arg == "smina":
            name = "Smina"
            cl = SminaDocking
        elif arg == "autodock4":
            name = "Autodock 4"
            cl = Autodock4Docking
    if name == None:
        nanome.util.Logs.error("Please pass the docking software to use as an argument: smina|autodock4")
        sys.exit(1)


    # Create the plugin, register Docking as the class to instantiate, and start listening
    plugin = nanome.Plugin(name + " Docking", "Run docking using " + name + ". Lets user choose the receptor, ligands, and diverse options", "Docking", True)
    plugin.set_plugin_class(cl)
    plugin.run('127.0.0.1', 8888)
