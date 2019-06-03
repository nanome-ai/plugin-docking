import nanome
from nanome.api.ui.image import Image as Image

class DockingMenu():
    def __init__(self, docking_plugin):
        self._plugin = docking_plugin
        self._selected_receptor = None
        self._selected_ligands = []
        self._selected_site = None
        self._exhaustiveness = 8
        self._modes = 9
        self._autobox_size = 4
        self._run_button = None
        self._align = True
        self._replace = False
        self._scoring_only = False
        self._tab = None
        self._autobox_enabled = True

    def _run_docking(self):
        if self._selected_receptor == None or len(self._selected_ligands) == 0:
            if self._autobox_enabled == True and self._selected_site == None:
                nanome.util.Logs.warning("Trying to run docking without having one receptor, one site and at least one ligand selected")
                return
        ligands = []
        for item in self._selected_ligands:
            ligands.append(item.complex)
        site = None
        if self._autobox_enabled:
            site = self._selected_site.complex
        self._plugin.run_docking(self._selected_receptor.complex, ligands, site)

    def disable_autobox(self):
        self._site_btn.unusable = True
        self._score_btn.unusable = True
        self._txt1.unusable = True # Doesn't do anything for now
        self._txt2.unusable = True
        self._txt3.unusable = True
        self._autobox_enabled = False
        self._plugin.update_menu(self._menu)

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
        self._receptor_checkmark.file_path = "checkmark.png"
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
            self._ligand_checkmark.file_path = "checkmark.png"
        else:
            self._ligand_checkmark.file_path = "none.png"
        self._plugin.update_content(self._ligand_checkmark)

    def site_pressed(self, button):
        lastSelected = self._selected_site
        if lastSelected != None:
            lastSelected.selected = False
            self._plugin.update_content(lastSelected)
        button.selected = True
        self._selected_site = button
        self._plugin.update_content(button)
        self._site_checkmark.file_path = "checkmark.png"
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
        self._receptor_checkmark.file_path = "none.png"
        self._ligand_checkmark.file_path = "none.png"
        self._site_checkmark.file_path = "none.png"

        for complex in complex_list:
            clone = self._complex_item_prefab.clone()
            ln_btn = clone.get_children()[0]
            btn = ln_btn.get_content()
            btn.set_all_text(complex.molecular.name)
            btn.complex = complex
            btn.register_pressed_callback(complex_pressed)
            self._complex_list.items.append(clone)

        self._plugin.update_menu(self._menu)

    def display_scoring_result(self, result):
        for molecule in result.molecules:
            clone = self._score_item_prefab.clone()
            ln_lbl = clone.get_children()[0]
            lbl = ln_lbl.get_content()
            lbl.text_value = molecule.molecular.name + " - " + molecule.molecular._associated["> <minimizedAffinity>"]
            self._score_list.items.append(clone)

        self._plugin.update_menu(self._menu)

    def build_menu(self):
        # defining callbacks
        def run_button_pressed_callback(button):
            if self._scoring_only:
                self._docking_param_panel.enabled = False
                self._score_panel.enabled = True
                self._score_list.items = []
                self._plugin.update_menu(self._menu)
            self._run_docking()

        def exhaustiveness_changed(input):
            try:
                self._exhaustiveness = int(input.input_text)
                nanome.util.Logs.debug("Exhaustiveness set to", self._exhaustiveness)
            except:
                self._exhaustiveness = 8
            if self._exhaustiveness <= 0:
                self._exhaustiveness = 8

        def modes_changed(input):
            try:
                self._modes = int(input.input_text)
                nanome.util.Logs.debug("Modes number set to", self._modes)
            except:
                self._modes = 9
            if self._modes <= 0:
                self._modes = 9

        def autobox_changed(input):
            try:
                self._autobox_size = int(input.input_text)
                nanome.util.Logs.debug("Autobox size set to", self._autobox_size)
            except:
                self._autobox_size = 4
            if self._autobox_size <= 0:
                self._autobox_size = 4

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
            elif button.text.value_idle == "Site":
                for item in self._complex_list.items:
                    btn = item.get_children()[0].get_content()
                    if btn == self._selected_site:
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

            self._plugin.update_menu(self._menu)

        def align_button_pressed_callback(button):
            self._align = not self._align
            button.selected = self._align
            self._plugin.update_content(button)

        def close_score_pressed_callback(button):
            self._docking_param_panel.enabled = True
            self._score_panel.enabled = False
            self._plugin.update_menu(self._menu)

        def scoring_button_pressed_callback(button):
            self._scoring_only = not self._scoring_only
            button.selected = self._scoring_only
            self._plugin.update_content(button)

        # Create a prefab that will be used to populate the lists
        self._complex_item_prefab = nanome.ui.LayoutNode()
        self._complex_item_prefab.layout_orientation = nanome.ui.LayoutNode.LayoutTypes.horizontal
        child = self._complex_item_prefab.create_child_node()
        child.forward_dist = 0.002
        child.add_new_button()

        self._score_item_prefab = nanome.ui.LayoutNode()
        self._score_item_prefab.layout_orientation = nanome.ui.LayoutNode.LayoutTypes.horizontal
        child = self._score_item_prefab.create_child_node()
        child.forward_dist = 0.002
        child.add_new_label()

        # loading menus
        menu = nanome.ui.Menu.io.from_json("_docking_menu.json")
        self._plugin.menu = menu

        # registering and saving special nodes

        # panels
        self._docking_param_panel = menu.root.find_node("LeftSide", True)
        self._score_panel = menu.root.find_node("LeftSideScore", True)

        # images
        self._receptor_checkmark = menu.root.find_node("ReceptorIcon", True).add_new_image("none.png")
        self._receptor_checkmark.scaling_option = Image.ScalingOptions.fit
        self._ligand_checkmark = menu.root.find_node("LigandIcon", True).add_new_image("none.png")
        self._ligand_checkmark.scaling_option = Image.ScalingOptions.fit
        self._site_checkmark = menu.root.find_node("SiteIcon", True).add_new_image("none.png")
        self._site_checkmark.scaling_option = Image.ScalingOptions.fit

        # texts
        self._txt1 = menu.root.find_node("ExhaustivenessInput", True).get_content()
        self._txt1.register_changed_callback(exhaustiveness_changed)

        self._txt2 = menu.root.find_node("ModesInput", True).get_content()
        self._txt2.register_changed_callback(modes_changed)

        self._txt3 = menu.root.find_node("AutoboxInput", True).get_content()
        self._txt3.register_changed_callback(autobox_changed)

        # buttons
        receptor_btn = menu.root.find_node("ReceptorButton", True).get_content()
        receptor_btn.register_pressed_callback(tab_button_pressed_callback)
        receptor_btn.selected = True
        self._tab = receptor_btn

        ligand_btn = menu.root.find_node("LigandButton", True).get_content()
        ligand_btn.register_pressed_callback(tab_button_pressed_callback)

        self._site_btn = menu.root.find_node("SiteButton", True).get_content()
        self._site_btn.register_pressed_callback(tab_button_pressed_callback)

        align_btn = menu.root.find_node("AlignButton", True).get_content()
        align_btn.register_pressed_callback(align_button_pressed_callback)
        align_btn.selected = True

        self._score_btn = menu.root.find_node("ScoringButton", True).get_content()
        self._score_btn.register_pressed_callback(scoring_button_pressed_callback)

        close_score_btn = menu.root.find_node("CloseScoreButton", True).get_content()
        close_score_btn.register_pressed_callback(close_score_pressed_callback)

        run_button = menu.root.find_node("RunButton", True).get_content()
        run_button.register_pressed_callback(run_button_pressed_callback)
        self._run_button = run_button

        # lists
        self._complex_list = menu.root.find_node("ComplexList", True).get_content()
        self._score_list = menu.root.find_node("ScoreList", True).get_content()
        self._receptor_tab = menu.root.find_node("ReceptorButton", True).get_content()
        self._ligand_tab = menu.root.find_node("LigandButton", True).get_content()
        self._site_tab = menu.root.find_node("SiteButton", True).get_content()

        # Update the menu
        self._menu = menu
        self._plugin.update_menu(menu)
        nanome.util.Logs.debug("Constructed plugin menu")
