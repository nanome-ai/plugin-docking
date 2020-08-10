import nanome
from nanome.api.ui.image import Image as Image
from nanome.util import Logs
import os
from nanome.api.ui import Dropdown,DropdownItem
from functools import partial


class DockingMenu():
    def __init__(self, docking_plugin):
        self._plugin = docking_plugin
        self._selected_receptor = None
        self._selected_ligands = []
        self._selected_site = None
        self._exhaustiveness = 10
        self._modes = 5
        self._autobox = 4
        self._run_button = None
        self._align = True
        self._replace = False
        self._scoring = False
        self._visual_scores = False
        self._tab = None
        self._autobox_enabled = True
    
    def get_params(self):
        params = {"exhaustiveness": None, "modes": None, "align": None, "replace": None, "scoring": None, "autobox": None}
        for key, value in params.items():
            newvalue = getattr(self, "_"+key)
            params[key] = newvalue
        return params

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
        params = {"exhaustiveness": None, "modes": None, "align": None, "replace": None, "scoring": None, "visual_scores": None, "autobox": None}
        for key, value in params.items():
            newvalue = getattr(self, "_"+key)
            params[key] = newvalue
        return params

    def _run_docking(self):
        if self._selected_receptor == None or len(self._selected_ligands) == 0:
            if self._autobox_enabled == True and self._selected_site == None:
                nanome.util.Logs.warning("Trying to run docking without having one receptor, one site and at least one ligand selected")
                return
        ligands = []
        for item in self._selected_ligands:
            # ligands.append(item.complex)
            ligands.append(item)
        site = None
        if self._autobox_enabled:
            # site = self._selected_site.complex 
            site = self._selected_site
        self.show_loading(True)
        self._plugin.run_docking(self._selected_receptor, ligands, site, self.get_params())

    def disable_autobox(self):
        self._site_btn.unusable = True
        self._score_btn.unusable = True
        #self._txt1.unusable = True # Doesn't do anything for now
        self._txt2.unusable = True
        #self._txt3.unusable = True
        self._autobox_enabled = False
        self._plugin.update_menu(self._menu)

    def make_plugin_usable(self, state=True):
        self._run_button.unusable = (not state) | self.refresh_run_btn_unusable(False)
        self._plugin.update_content(self._run_button)
    
    def show_loading(self, show = False):
        if show:
            self.ln_run_button.enabled = False
            self.ln_loading_bar.enabled = True
            Logs.debug("loading bar show=========================")
        else:
            self.ln_run_button.enabled = True
            self.ln_loading_bar.enabled = False
        #self._plugin.update_node(self.ln_loading_bar)
        self._plugin.update_menu(self._menu)

    def receptor_pressed(self, button):
        lastSelected = self._selected_receptor
        if lastSelected != None:
            lastSelected.selected = False
            self._plugin.update_content(lastSelected)
        button.selected = True
        self._selected_receptor = button
        self._plugin.update_content(button)
        #self._receptor_checkmark.file_path = os.path.join(os.path.dirname(__file__), 'checkmark.png')
        self._plugin.update_content(self._receptor_checkmark)

        self.refresh_run_btn_unusable()

    def ligand_pressed(self, button):
        if button.selected == False:
            button.selected = True
            self._selected_ligands.append(button)
        else:
            button.selected = False
            self._selected_ligands.remove(button)
        self._plugin.update_content(button)
        if len(self._selected_ligands) != 0:
            self._ligand_checkmark.file_path = os.path.join(os.path.dirname(__file__), 'checkmark.png')
        else:
            self._ligand_checkmark.file_path = os.path.join(os.path.dirname(__file__), 'none.png')
        self._plugin.update_content(self._ligand_checkmark)

        self.refresh_run_btn_unusable()

    def site_pressed(self, button):
        lastSelected = self._selected_site
        if lastSelected != None:
            lastSelected.selected = False
            self._plugin.update_content(lastSelected)
        button.selected = True
        self._selected_site = button
        self._plugin.update_content(button)
        # self._site_checkmark.file_path = os.path.join(os.path.dirname(__file__), 'checkmark.png')
        # self._plugin.update_content(self._site_checkmark)

        self.refresh_run_btn_unusable()

    def refresh_run_btn_unusable(self, update=True):
        site_requirement_met = self._selected_site != None or not self._plugin._calculations.requires_site
        if self._selected_receptor != None and len(self._selected_ligands) > 0 and site_requirement_met:
            self._run_button.text.value_unusable = "Running..."
            self._run_button.unusable = False
        else:
            self._run_button.text.value_unusable = "Run"
            self._run_button.unusable = True
        if update:
            self._plugin.update_content(self._run_button)

        return self._run_button.unusable

    def change_complex_list(self, complex_list):
        def complex_pressed(button):
            if self._tab.text.value_idle == "Receptor":
                self.receptor_pressed(button)
            elif self._tab.text.value_idle == "Ligand":
                self.ligand_pressed(button)
            elif self._tab.text.value_idle == "Site":
                self.site_pressed(button)

        self.reset(update_menu=False)

        # for complex in complex_list:
        #     clone = self._complex_item_prefab.clone()
        #     ln_btn = clone.get_children()[0]
        #     btn = ln_btn.get_content()
        #     btn.set_all_text(complex.full_name)
        #     btn.complex = complex
        #     btn.register_pressed_callback(complex_pressed)
        #     self._complex_list.items.append(clone)
        ligand_list = []
        receptor_list = []
        site_list = []
        for complex in complex_list:
            dd_item1 = DropdownItem()
            dd_item2 = DropdownItem()
            dd_item3 = DropdownItem()
            dd_item1.complex = complex
            dd_item2.complex = complex
            dd_item3.complex = complex
            dd_item1._name = complex.full_name
            dd_item2._name = complex.full_name
            dd_item3._name = complex.full_name
            ligand_list.append(dd_item1)
            receptor_list.append(dd_item2)
            site_list.append(dd_item3)

        self._ligand_dropdown.items = ligand_list
        self._receptor_dropdown.items = receptor_list
        self._site_dropdown.items = site_list
        self._ligand_dropdown.register_item_clicked_callback(partial(self.handle_dropdown_pressed,self._selected_ligands,'ligand'))
        self._receptor_dropdown.register_item_clicked_callback(partial(self.handle_dropdown_pressed,self._selected_receptor,'receptor'))
        self._site_dropdown.register_item_clicked_callback(partial(self.handle_dropdown_pressed,self._selected_site,'site'))
        self._plugin.update_menu(self._menu)

    def display_scoring_result(self, result):
        self.reset()

        for molecule in result.molecules:
            clone = self._score_item_prefab.clone()
            ln_lbl = clone.get_children()[0]
            lbl = ln_lbl.get_content()
            lbl.text_value = molecule.name + " - " + molecule._associated["> <minimizedAffinity>"]
            self._score_list.items.append(clone)

    def reset(self, update_menu=True):
        self._selected_receptor = None
        self._selected_ligands = []
        self._selected_site = None
        self._complex_list.items = []
        # self._receptor_checkmark.file_path = os.path.join(os.path.dirname(__file__), 'none.png')
        # self._ligand_checkmark.file_path = os.path.join(os.path.dirname(__file__), 'none.png')
        # self._site_checkmark.file_path = os.path.join(os.path.dirname(__file__), 'none.png')

        self.make_plugin_usable()
        self._plugin.update_menu(self._menu)

    def open_menu(self):
        self._plugin.menu = self._menu
        self._plugin.menu.enabled = True
        self._plugin.update_menu(self._plugin.menu)

    def handle_dropdown_pressed(self,docking_component,component_name,dropdown,item):
        if component_name == 'ligand':
            if item.complex.index not in [x.index for x in self._selected_ligands]:
                self._selected_ligands.append(item.complex)
            if len(self._selected_ligands) > 1:
                self._ligand_txt._text_value = 'Multiple'
            else:
                self._ligand_txt._text_value = item.complex.full_name if len(item.complex.full_name) <= 4 else item.complex.full_name[:8]+'...'
        elif component_name == 'receptor':
            self._selected_receptor = item.complex
            self._receptor_txt._text_value = item.complex.full_name if len(item.complex.full_name) <= 4 else item.complex.full_name[:8]+'...'
        elif component_name == 'site':
            self._selected_site = item.complex
            
        Logs.debug("dropdown is: ",dropdown)
        Logs.debug("item is: ",item)
        Logs.debug("docking component is ",docking_component)
        Logs.debug("selected ligand is ",self._selected_ligands)
        Logs.debug("selected receptor is ",self._selected_receptor)
        Logs.debug("selected site is: ",self._selected_site)
        Logs.debug()
        self.update_icons()
        self.refresh_run_btn_unusable()
        self._plugin.update_menu(self._menu)
 
    def update_icons(self):
        if self._selected_ligands:
            self._ligand_icon._file_path = self.ligand_white_path
        else:
            self._ligand_icon._file_path = self.ligand_gray_path


        if self._selected_receptor:
            self._receptor_icon._file_path = self.receptor_white_path
        else:
            self._receptor_icon._file_path = self.receptor_gray_path

    def build_menu(self):
        # defining callbacks
        def run_button_pressed_callback(button):
            if self._scoring:
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
                self._autobox = int(input.input_text)
                nanome.util.Logs.debug("Autobox size set to", self._autobox)
            except:
                self._autobox = 4
            if self._autobox <= 0:
                self._autobox = 4

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
            self._scoring = not self._scoring
            button.selected = self._scoring
            self._plugin.update_content(button)

        def visual_scores_button_pressed_callback(button):
            self._visual_scores = not self._visual_scores
            button.selected = self._visual_scores
            self._plugin.update_content(button)



        # Create a prefab that will be used to populate the lists
        self._complex_item_prefab = nanome.ui.LayoutNode()
        self._complex_item_prefab.layout_orientation = nanome.ui.LayoutNode.LayoutTypes.horizontal
        child = self._complex_item_prefab.create_child_node()
        child.add_new_button()

        self._score_item_prefab = nanome.ui.LayoutNode()
        self._score_item_prefab.layout_orientation = nanome.ui.LayoutNode.LayoutTypes.horizontal
        child = self._score_item_prefab.create_child_node()
        # child.forward_dist = 0.002
        child.add_new_label()

        # loading menus
        menu = nanome.ui.Menu.io.from_json(os.path.join(os.path.dirname(__file__), '_docking_menu_new.json'))
        self._plugin.menu = menu

        # registering and saving special nodes

        # panels
        self._docking_param_panel = menu.root.find_node("LeftSide", True)
        self._score_panel = menu.root.find_node("LeftSideScore", True)
        self._panel_separator = menu.root.find_node("MiddleLine",True)
        

        # images
        self.can_dock_path = os.path.join(os.path.dirname(__file__), 'can_dock.png')
        self.cannot_dock_path = os.path.join(os.path.dirname(__file__), 'cannot_dock.png')
        self._check_arrow = menu.root.find_node("CheckArrow",True).add_new_image(self.cannot_dock_path)

        self.ligand_gray_path = os.path.join(os.path.dirname(__file__), 'ligand_gray.png')
        self.ligand_white_path = os.path.join(os.path.dirname(__file__), 'ligand_white.png')
        self._ligand_icon = menu.root.find_node("LigandIcon",True).add_new_image(self.ligand_gray_path)

        self.receptor_gray_path = os.path.join(os.path.dirname(__file__), 'receptor_gray.png')
        self.receptor_white_path = os.path.join(os.path.dirname(__file__), 'receptor_white.png')
        self._receptor_icon = menu.root.find_node("ReceptorIcon",True).add_new_image(self.receptor_gray_path)

        # none_path = os.path.join(os.path.dirname(__file__), 'none.png')
        # self._receptor_checkmark = menu.root.find_node("ReceptorIcon", True).add_new_image(none_path)
        # self._receptor_checkmark.scaling_option = Image.ScalingOptions.fit
        # self._ligand_checkmark = menu.root.find_node("LigandIcon", True).add_new_image(none_path)
        # self._ligand_checkmark.scaling_option = Image.ScalingOptions.fit
        # self._site_checkmark = menu.root.find_node("SiteIcon", True).add_ new_image(none_path)
        # self._site_checkmark.scaling_option = Image.ScalingOptions.fit

        # texts
        #self._txt1 = menu.root.find_node("ExhaustivenessInput", True).get_content()
        #self._txt1.register_changed_callback(exhaustiveness_changed)

        self._txt2 = menu.root.find_node("ModesInput", True).get_content()
        self._txt2.register_changed_callback(modes_changed)
        self._ligand_txt = menu.root.find_node("LigandName").get_content()
        self._receptor_txt = menu.root.find_node("ReceptorName").get_content()
        #self._txt3 = menu.root.find_node("AutoboxInput", True).get_content()
        #self._txt3.register_changed_callback(autobox_changed)

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

        self._score_btn = menu.root.find_node("ScoringButton").get_content()
        self._score_btn.register_pressed_callback(scoring_button_pressed_callback)

        self._display_score_btn = menu.root.find_node("VisualScoresButton").get_content()
        self._display_score_btn.register_pressed_callback(visual_scores_button_pressed_callback)

        close_score_btn = menu.root.find_node("CloseScoreButton").get_content()
        close_score_btn.register_pressed_callback(close_score_pressed_callback)

        self.ln_run_button = menu.root.find_node("RunButton")
        run_button = self.ln_run_button.get_content()
        run_button.register_pressed_callback(run_button_pressed_callback)
        self._run_button = run_button
        self._run_button.enabled = False
        self.refresh_run_btn_unusable()

        # loading bar
        self.ln_loading_bar = menu.root.find_node("LoadingBar")
        self.loading_bar = self.ln_loading_bar.get_content()
        self.loading_bar.description = "      Loading...          "

        # lists
        self._complex_list = menu.root.find_node("ComplexList").get_content()
        self._score_list = menu.root.find_node("ScoreList").get_content()
        self._receptor_tab = menu.root.find_node("ReceptorButton").get_content()
        self._ligand_tab = menu.root.find_node("LigandButton").get_content()
        self._site_tab = menu.root.find_node("SiteButton").get_content()

        # dropdown
        self._ligand_dropdown = menu.root.find_node("LigandDropdown").add_new_dropdown()
        self._receptor_dropdown = menu.root.find_node("ReceptorDropdown").add_new_dropdown()
        self._site_dropdown = menu._root.find_node("SiteDropdown").add_new_dropdown()

        # Update the menu
        self._menu = menu
        self._plugin.update_menu(menu)
        nanome.util.Logs.debug("Constructed plugin menu")
