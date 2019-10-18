import nanome
from ._DockingCalculations import DockingCalculations as Smina
from ._DockingCalculationsAutodock4 import DockingCalculations as Autodock4
from ._DockingMenu import DockingMenu
import sys

__metaclass__ = type


class Docking(nanome.PluginInstance):
    def __init__(self):
        self._menu = DockingMenu(self)
        self._calculations = None
        self._autobox = True

    # Function called when Nanome connects to the Plugin, after its instantiation
    def start(self):
        self._menu.build_menu()
        if self._autobox == False:
            self._menu.disable_autobox()
        self.request_complex_list(self.on_complex_list_received)

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

    def on_complex_added(self):
        self.request_complex_list(self.on_complex_list_received)

    def on_complex_removed(self):
        self.request_complex_list(self.on_complex_list_received)

    def open_menu(self):
        menu = self.menu
        menu.enabled = True
        self.update_menu(menu)

    def make_plugin_usable(self):
        self._menu.make_plugin_usable()

    # Function called when Nanome returns the complex list after a request
    def on_complex_list_received(self, complexes):
        self._menu.change_complex_list(complexes)

    def run_docking(self, receptor, ligand_list, site):
        has_site = site != None

        def on_complexes_received(complexes):
            receptor = complexes[0]
            self._receptor = receptor
            Docking.convert_atoms_to_absolute_position(receptor)
            starting_lig_idx = 1
            site = None
            if has_site:
                site = complexes[1]
                self._site = site
                self._site_wtc_mat = site.get_workspace_to_complex_matrix()
                Docking.convert_atoms_to_absolute_position(site)
                starting_lig_idx = 2
            ligands = nanome.structure.Complex()
            for ligand in complexes[starting_lig_idx:]:
                Docking.convert_atoms_to_absolute_position(ligand)
                for molecule in ligand.molecules:
                    ligands.add_molecule(molecule)
            self._calculations.start_docking(receptor, ligands, site, **params)

        request_list = [receptor.index]
        if has_site:
            request_list.append(site.index)
        request_list += [x.index for x in ligand_list]
        self.request_complexes(request_list, on_complexes_received)

    @staticmethod
    def convert_atoms_to_absolute_position(complex):
        mat = complex.get_complex_to_workspace_matrix()
        for atom in complex.atoms:
            atom.position = mat * atom.position

    @staticmethod
    def convert_atoms_to_relative_position(complex, mat):
        # mat = reference.get_workspace_to_complex_matrix()
        for atom in complex.atoms:
            atom.position = mat * atom.position

    # Function called every update tick of the Plugin
    def update(self):
        self._calculations.update()

    def add_result_to_workspace(self, results, align=False):
        for complex in results:
            # reference = site if site != None else self._receptor
            # print("reference:", [chain.name for chain in reference.chains])
            Docking.convert_atoms_to_relative_position(complex, self._site_wtc_mat)

            if align:
                # get global atom pos
                mat = complex.get_complex_to_workspace_matrix()
                global_pos = [mat * atom.position for atom in complex.atoms]
                # align bounding box
                complex.position = self._site.position
                complex.rotation = self._site.rotation
                # restore atom pos
                mat = complex.get_workspace_to_complex_matrix()
                for (atom, pos) in zip(complex.atoms, global_pos):
                    atom.position = mat * pos

                complex.boxed = True

        self.add_to_workspace(results)

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
        self._autobox = False


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
    if name == None:
        nanome.util.Logs.error("Please pass the docking software to use as an argument: smina|autodock4")
        sys.exit(1)

    # Create the plugin, register Docking as the class to instantiate, and start listening
    plugin = nanome.Plugin(name + " Docking", "Run docking using " + name + ". Lets user choose the receptor, ligands, and diverse options", "Docking", True)
    plugin.set_plugin_class(cl)
    plugin.run('127.0.0.1', 8888)


if __name__ == "__main__":
    main()
