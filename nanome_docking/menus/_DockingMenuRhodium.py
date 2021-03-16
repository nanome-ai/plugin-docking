import nanome
from nanome.api.ui.image import Image as Image

import os

BASE_DIR = os.path.dirname(__file__)
ICONS_DIR = os.path.join(BASE_DIR, 'icons')
ICONS = { icon.rsplit('.')[0]: os.path.join(ICONS_DIR, icon) for icon in os.listdir(ICONS_DIR) }

class DockingMenuRhodium():
    def __init__(self, docking_plugin):
        self._plugin = docking_plugin
        self._selected_receptor = None
        self._selected_ligands = []
        self._selected_site = None
        self._grid_resolution = 2.5
        self._poses = 128
        self._rotamer = 6
        self._run_button = None
        self._align = True
        self._ignore_hetatoms = False
        self._tab = None

    def is_ready_for_docking(self):
        return self._selected_receptor != None and len(self._selected_ligands) > 0

    def get_receptor(self):
        return self._selected_receptor.complex

    def get_ligands(self):
        ligands = []
        for item in self._selected_ligands:
            ligands.append(item.complex)
        return ligands

    def get_site(self):
        if self._selected_site == None:
            return None
        return self._selected_site.complex

    def get_params(self):
        return { "receptor":self._selected_receptor,
                "ligands":self._selected_ligands,
                "site":self._selected_site,
                "grid_resolution":self._grid_resolution,
                "poses":self._poses,
                "rotamers":self._rotamer,
                "align":self._align,
                "ignore_hetatoms":self._ignore_hetatoms
            }

    def _run_docking(self):
        if self._selected_receptor == None or len(self._selected_ligands) == 0:
            nanome.util.Logs.warning("Trying to run docking without having one receptor, and at least one ligand selected")
            self._plugin.send_notification(nanome.util.enums.NotificationTypes.error, "Please select a receptor, and at least one ligand")
            return
        ligands = []
        for item in self._selected_ligands:
            ligands.append(item.complex)
        self._plugin.run_docking(self.get_receptor(), self.get_ligands(), self.get_site(), self.get_params())

    def make_plugin_usable(self, state=True):
        self._run_button.unusable = not state
        self._plugin.update_content(self._run_button)

    def receptor_pressed(self, button):
        lastSelected = self._selected_receptor
        if lastSelected != None:
            lastSelected.selected = False
            self._plugin.update_content(lastSelected)
        button.selected = True
        self._selected_receptor = button
        self._plugin.update_content(button)
        self._receptor_checkmark.file_path = ICONS['checkmark']
        self._plugin.update_content(self._receptor_checkmark)

    def ligand_pressed(self, button):
        if button.selected == False:
            button.selected = True
            self._selected_ligands.append(button)
        else:
            button.selected = False
            self._selected_ligands.remove(button)
        self._plugin.update_content(button)
        if len(self._selected_ligands) != 0:
            self._ligand_checkmark.file_path = ICONS['checkmark']
        else:
            self._ligand_checkmark.file_path = ICONS['none']
        self._plugin.update_content(self._ligand_checkmark)

    def site_pressed(self, button):
        lastSelected = self._selected_site
        if lastSelected == button:
            button.selected = False
            self._plugin.update_content(button)
            self._selected_site = None
            self._site_checkmark.file_path = ICONS['none']
            self._plugin.update_content(self._site_checkmark)
            return
        elif lastSelected != None:
            lastSelected.selected = False
            self._plugin.update_content(lastSelected)
        button.selected = True
        self._selected_site = button
        self._plugin.update_content(button)
        self._site_checkmark.file_path = ICONS['checkmark']
        self._plugin.update_content(self._site_checkmark)

    def change_complex_list(self, complex_list):
        def complex_pressed(button):
            if self._tab.text.value_idle == "Receptor":
                self.receptor_pressed(button)
            elif self._tab.text.value_idle == "Ligand":
                self.ligand_pressed(button)
            elif self._tab.text.value_idle == "Site":
                self.site_pressed(button)

        self._selected_receptor = None
        self._selected_ligands = []
        self._selected_site = None
        self._complex_list.items = []
        self._receptor_checkmark.file_path = ICONS['none']
        self._ligand_checkmark.file_path = ICONS['none']
        self._site_checkmark.file_path = ICONS['none']

        for complex in complex_list:
            clone = self._complex_item_prefab.clone()
            ln_btn = clone.get_children()[0]
            btn = ln_btn.get_content()
            btn.set_all_text(complex.name)
            btn.complex = complex
            btn.register_pressed_callback(complex_pressed)
            self._complex_list.items.append(clone)

        self._plugin.update_menu(self._menu)

    def build_menu(self):
        # defining callbacks
        def run_button_pressed_callback(button):
            self._run_docking()

        def grid_resolution_changed(input):
            try:
                self._grid_resolution = float(input.input_text)
                nanome.util.Logs.debug("Grid resolution set to", self._grid_resolution)
            except:
                self._grid_resolution = 1.4

        def poses_changed(input):
            try:
                self._poses = int(input.input_text)
                nanome.util.Logs.debug("Poses count set to", self._poses)
            except:
                self._poses = 96
            if self._poses <= 0:
                self._poses = 96

        def rotamer_changed(input):
            try:
                self._rotamer = int(input.input_text)
                nanome.util.Logs.debug("Rotamer count set to", self._rotamer)
            except:
                self._rotamer = 3
            if self._rotamer <= 0:
                self._rotamer = 3

        def tab_button_pressed_callback(button):
            if self._tab == button:
                return

            self._tab.selected = False
            button.selected = True
            self._tab = button

            if button.text.value_idle == "Receptor":
                for item in self._complex_list.items:
                    btn = item.get_children()[0].get_content()
                    if btn == self._selected_receptor:
                        btn.selected = True
                    else:
                        btn.selected = False
            elif button.text.value_idle == "Ligand":
                for item in self._complex_list.items:
                    btn = item.get_children()[0].get_content()
                    if btn in self._selected_ligands:
                        btn.selected = True
                    else:
                        btn.selected = False
            elif button.text.value_idle == "Site":
                for item in self._complex_list.items:
                    btn = item.get_children()[0].get_content()
                    if btn == self._selected_site:
                        btn.selected = True
                    else:
                        btn.selected = False

            self._plugin.update_menu(self._menu)

        def align_button_pressed_callback(button):
            self._align = not self._align
            button.selected = self._align
            self._plugin.update_content(button)

        def ignore_hetatoms_pressed_callback(button):
            self._ignore_hetatoms = not self._ignore_hetatoms
            button.selected = self._ignore_hetatoms
            self._plugin.update_content(button)

        # Create a prefab that will be used to populate the lists
        self._complex_item_prefab = nanome.ui.LayoutNode()
        self._complex_item_prefab.layout_orientation = nanome.ui.LayoutNode.LayoutTypes.horizontal
        child = self._complex_item_prefab.create_child_node()
        child.forward_dist = 0.002
        child.add_new_button()

        # loading menus
        menu = nanome.ui.Menu.io.from_json(os.path.join(BASE_DIR, 'jsons', '_docking_menu_rhodium.json'))
        self._plugin.menu = menu

        # images
        self._receptor_checkmark = menu.root.find_node("ReceptorIcon").add_new_image(ICONS['none'])
        self._receptor_checkmark.scaling_option = Image.ScalingOptions.fit
        self._ligand_checkmark = menu.root.find_node("LigandIcon").add_new_image(ICONS['none'])
        self._ligand_checkmark.scaling_option = Image.ScalingOptions.fit
        self._site_checkmark = menu.root.find_node("SiteIcon").add_new_image(ICONS['none'])
        self._site_checkmark.scaling_option = Image.ScalingOptions.fit
        logo_image = menu.root.find_node("LogoImage").add_new_image(icons['swri_logo'])
        logo_image.scaling_option = Image.ScalingOptions.fit

        # texts
        txt1 = menu.root.find_node("GridResolutionInput").get_content()
        txt1.register_changed_callback(grid_resolution_changed)

        txt2 = menu.root.find_node("PosesCountInput").get_content()
        txt2.register_changed_callback(poses_changed)

        txt3 = menu.root.find_node("RotamerCountInput").get_content()
        txt3.register_changed_callback(rotamer_changed)

        # buttons
        receptor_btn = menu.root.find_node("ReceptorButton").get_content()
        receptor_btn.register_pressed_callback(tab_button_pressed_callback)
        receptor_btn.selected = True
        self._tab = receptor_btn

        ligand_btn = menu.root.find_node("LigandButton").get_content()
        ligand_btn.register_pressed_callback(tab_button_pressed_callback)

        site_btn = menu.root.find_node("SiteButton").get_content()
        site_btn.register_pressed_callback(tab_button_pressed_callback)

        align_btn = menu.root.find_node("AlignButton").get_content()
        align_btn.register_pressed_callback(align_button_pressed_callback)
        align_btn.selected = True

        run_button = menu.root.find_node("RunButton").get_content()
        run_button.register_pressed_callback(run_button_pressed_callback)
        self._run_button = run_button

        # lists
        self._complex_list = menu.root.find_node("ComplexList").get_content()
        self._receptor_tab = menu.root.find_node("ReceptorButton").get_content()
        self._ligand_tab = menu.root.find_node("LigandButton").get_content()
        self._site_tab = menu.root.find_node("SiteButton").get_content()

        # Update the menu
        self._menu = menu
        self._plugin.update_menu(menu)
        nanome.util.Logs.debug("Constructed plugin menu")
