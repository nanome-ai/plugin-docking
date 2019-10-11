import nanome
from ._DockingCalculations import DockingCalculations as Smina
from ._DockingCalculationsAutodock4 import DockingCalculations as Autodock4
from ._DockingCalculationsRhodium import DockingCalculations as Rhodium
from ._DockingMenu import DockingMenu
from ._DockingMenuRhodium import DockingMenuRhodium
import sys

__metaclass__ = type


class Docking(nanome.PluginInstance):
    def __init__(self):
        self._menu = None
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
        menu = self._menu
        # If menu doesn't have Receptor and Ligands selected, open it
        # Else, just start docking
        if menu.is_ready_for_docking() == False:
            self.open_menu()
        else:
            self.run_docking(menu.get_receptor(), menu.get_ligands(), menu.get_site(), menu.get_params())

    # Called when user click on the "Advanced Settings" button in Nanome
    def on_advanced_settings(self):
        nanome.util.Logs.debug("Advanced Settings")
        self.open_menu()

    # Called when a complex is added to the workspace in Nanome
    def on_complex_added(self):
        self.request_complex_list(self.on_complex_list_received)

    # Called when a complex is removed from the workspace in Nanome
    def on_complex_removed(self):
        self.request_complex_list(self.on_complex_list_received)

    def open_menu(self):
        menu = self.menu
        menu.enabled = True
        self.update_menu(menu)

    def make_plugin_usable(self):
        self._menu.make_plugin_usable()

    def on_complex_list_received(self, complexes):
        self._menu.change_complex_list(complexes)

    def run_docking(self, receptor, ligands, site, params):
        # Change the plugin to be "unusable"
        if self._menu._run_button.unusable == True:
            return
        self._menu.make_plugin_usable(False)

        has_site = site != None

        def on_complexes_received(complexes):
            # When deep complexes data are received, unpack them and prepare ligand for docking
            receptor = complexes[0]
            self._receptor = receptor
            Docking.convert_atoms_to_absolute_position(receptor)
            starting_lig_idx = 1
            site = None
            if has_site:
                site = complexes[1]
                Docking.convert_atoms_to_absolute_position(site)
                starting_lig_idx = 2
            ligands = nanome.structure.Complex()
            for ligand in complexes[starting_lig_idx:]:
                Docking.convert_atoms_to_absolute_position(ligand)
                for molecule in ligand.molecules:
                    ligands.add_molecule(molecule)
            self._calculations.start_docking(receptor, ligands, site, **params)

        # Request complexes to Nanome in this order: [receptor, site (if any), ligand, ligand,...]
        request_list = [receptor.index]
        if has_site:
            request_list.append(site.index)
        request_list += [x.index for x in ligands]
        self.request_complexes(request_list, on_complexes_received)

    @staticmethod
    def convert_atoms_to_absolute_position(complex):
        mat = complex.transform.get_complex_to_workspace_matrix()
        for atom in complex.atoms:
            atom.molecular.position = mat * atom.molecular.position

    @staticmethod
    def convert_atoms_to_relative_position(complex, reference):
        mat = reference.transform.get_workspace_to_complex_matrix()
        for atom in complex.atoms:
            atom.molecular.position = mat * atom.molecular.position

    # Called every update tick of the Plugin
    def update(self):
        self._calculations.update()

    def add_result_to_workspace(self, results):
        for complex in results:
            Docking.convert_atoms_to_relative_position(complex, self._receptor)
        self.add_to_workspace(results)

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
        nanome.util.Logs.error("Please pass the docking software to use as an argument: smina|autodock4|rhodium")
        sys.exit(1)

    # Create the plugin, register Docking as the class to instantiate, and start listening
    plugin = nanome.Plugin(name + " Docking", "Run docking using " + name + ". Lets user choose the receptor, ligands, and diverse options", "Docking", True)
    plugin.set_plugin_class(cl)
    plugin.run('127.0.0.1', 8888)


if __name__ == "__main__":
    main()
