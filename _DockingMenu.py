import nanome

class DockingMenu():
    def __init__(self, docking_plugin):
        self._plugin = docking_plugin
        self._selected_receptor = None # button with a complex
        self._selected_ligands = [] # List of complex
        self._selected_site = None # button with a complex
        self._complex_list = []
        self._exhaustiveness = 8
        self._modes = 9
        self._autobox_size = 4
        self._run_button = None
        self._align = False
        self._replace = False
        self._scoring_only = False

    def _request_refresh(self):
        self._plugin.request_refresh()

    def _run_docking(self):
        if self._selected_receptor == None or self._selected_site == None or self._selected_ligands.count == 0:
            nanome.util.Logs.warning("Trying to run docking without having one receptor, one site and at least one ligand selected")
            return
        self._plugin.run_docking(self._selected_receptor.complex, self._selected_ligands, self._selected_site.complex)

    def _find_complex(self, name):
        for complex in self._complex_list:
            if complex.molecular.name == name:
                return complex
        return None

    def make_plugin_usable(self, state = True):
        self._run_button.unusable = not state
        self._plugin.update_content(self._run_button)

    def change_complex_list(self, complex_list):
        def receptor_pressed(button):
            complex = button.complex
            if complex == None:
                nanome.util.Logs.error("Couldn't retrieve a complex from its button text")
                return
            lastSelected = self._selected_receptor
            if lastSelected != None:
                lastSelected.selected = False
                self._plugin.update_content(lastSelected)
            button.selected = True
            self._selected_receptor = button
            self._plugin.update_content(button)
            
        def ligand_pressed(button):
            complex = button.complex
            if complex == None:
                nanome.util.Logs.error("Couldn't retrieve a complex from its button text")
                return

            if button.selected == False:
                button.selected = True
                self._selected_ligands.append(complex)
            else:
                button.selected = False
                self._selected_ligands.remove(complex)
            self._plugin.update_content(button)

        def site_pressed(button):
            complex = button.complex
            if complex == None:
                nanome.util.Logs.error("Couldn't retrieve a complex from its button text")
                return
            lastSelected = self._selected_site
            if lastSelected != None:
                lastSelected.selected = False
                self._plugin.update_content(lastSelected)
            button.selected = True
            self._selected_site = button
            self._plugin.update_content(button)

        self._complex_list = complex_list
        self._selected_receptor = None
        self._selected_ligands = []
        self._selected_site = None
        self._receptor_list.items = []
        self._ligands_list.items = []
        self._site_list.items = []

        for complex in complex_list:
            nanome.util.Logs.debug("Adding to list:", complex.molecular.name)
            clone = self._complex_item_prefab.clone()
            ln_btn = clone.get_children()[0]
            btn = ln_btn.get_content()
            btn.set_all_text(complex.molecular.name)
            btn.complex = complex
            btn.register_pressed_callback(receptor_pressed)
            self._receptor_list.items.append(clone)
            clone1 = clone.clone()
            ln_btn = clone1.get_children()[0]
            btn = ln_btn.get_content()
            btn.complex = complex
            btn.register_pressed_callback(ligand_pressed)
            self._ligands_list.items.append(clone1)
            clone2 = clone.clone()
            ln_btn = clone2.get_children()[0]
            btn = ln_btn.get_content()
            btn.complex = complex
            btn.register_pressed_callback(site_pressed)
            self._site_list.items.append(clone2)

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
        #defining callbacks
        def refresh_button_pressed_callback(button):
            self._request_refresh()

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
        child.name = "button_node"
        child.add_new_button()

        self._score_item_prefab = nanome.ui.LayoutNode()
        self._score_item_prefab.layout_orientation = nanome.ui.LayoutNode.LayoutTypes.horizontal
        child = self._score_item_prefab.create_child_node()
        child.name = "label_node"
        child.add_new_label()

        # loading menus
        menu = nanome.ui.Menu.io.from_json("_docking_menu.json")
        nanome.ui.Menu.set_plugin_menu(menu)

        # registering and saving special nodes

        # panels
        self._docking_param_panel = menu.root.find_node("docking_panel", True)
        self._score_panel = menu.root.find_node("score_panel", True)

        # texts
        txt1 = menu.root.find_node("txt1", True).get_content()
        txt1.register_changed_callback(exhaustiveness_changed)

        txt2 = menu.root.find_node("txt2", True).get_content()
        txt2.register_changed_callback(modes_changed)

        txt3 = menu.root.find_node("txt3", True).get_content()
        txt2.register_changed_callback(autobox_changed)

        align_btn = menu.root.find_node("Align", True).get_content()
        align_btn.register_pressed_callback(align_button_pressed_callback)

        score_btn = menu.root.find_node("Score", True).get_content()
        score_btn.register_pressed_callback(scoring_button_pressed_callback)

        close_score_btn = menu.root.find_node("Close_Score", True).get_content()
        close_score_btn.register_pressed_callback(close_score_pressed_callback)

        refresh_button = menu.root.find_node("Refresh", True).get_content()
        refresh_button.register_pressed_callback(refresh_button_pressed_callback)

        run_button = menu.root.find_node("Run", True).get_content()
        run_button.register_pressed_callback(run_button_pressed_callback)
        self._run_button = run_button

        # lists
        self._receptor_list = menu.root.find_node("left_list", True).get_content()
        self._ligands_list = menu.root.find_node("center_list", True).get_content()
        self._site_list = menu.root.find_node("right_list", True).get_content()
        self._score_list = menu.root.find_node("score_list", True).get_content()

        # Update the menu
        self._menu = menu
        self._plugin.update_menu(self._menu)
        self._request_refresh()
        nanome.util.Logs.debug("Constructed plugin menu")