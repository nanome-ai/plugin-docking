import asyncio
import os
from functools import partial

import nanome
from nanome.util import Logs, async_callback, Vector3
from nanome.api.ui import DropdownItem
from nanome.api.shapes import Sphere, Shape

BASE_DIR = os.path.dirname(__file__)
ICONS_DIR = os.path.join(BASE_DIR, 'icons')
ICONS = {icon.rsplit('.')[0]: os.path.join(ICONS_DIR, icon) for icon in os.listdir(ICONS_DIR)}


class DockingMenu():

    def __init__(self, docking_plugin):
        self._plugin = docking_plugin
        self._selected_receptor = None
        self._selected_ligands = []
        self._selected_site = None
        self._modes = 5
        self._autobox = 4
        self._run_button = None
        self._align = True
        self._replace = False
        self._scoring = False
        self._visual_scores = False
        self._tab = None

        # loading menus
        self._menu = nanome.ui.Menu.io.from_json(os.path.join(BASE_DIR, 'jsons', '_docking_menu.json'))

        # Run button
        self.ln_run_button = self._menu.root.find_node("RunButton")
        self._run_button = self.ln_run_button.get_content()
        # loading bar
        self.ln_loading_bar = self._menu.root.find_node("LoadingBar")
        self.loading_bar = self.ln_loading_bar.get_content()
        self.loading_bar.description = "    Loading...          "

        algo_name = self._plugin.__class__.__name__.split('Docking')[0]
        self._menu.title = f'{algo_name} Docking'

    def get_receptor(self):
        return self._selected_receptor.complex

    def get_ligands(self):
        ligands = []
        for item in self._selected_ligands:
            ligands.append(item.complex)
        return ligands

    def get_site(self):
        if self._selected_site is None:
            return None
        return self._selected_site.complex

    def get_params(self):
        """Collect parameters from this menu and the Settings Menu."""
        params = {
            "exhaustiveness": None,
            "modes": None,
            "align": None,
            "replace": None,
            "scoring": None,
            "visual_scores": None,
            "autobox": None
        }
        settings_menu = self._plugin.settings_menu
        for key in params.keys():
            attr_key = f'_{key}'
            if hasattr(self, attr_key):
                newvalue = getattr(self, attr_key)
            else:
                # Get value from Settings menu
                newvalue = getattr(settings_menu, attr_key)
            params[key] = newvalue
        return params

    def create_complex_dropdown_items(self, complex_list):
        ddi_list = []
        for comp in complex_list:
            ddi = DropdownItem(comp.full_name)
            ddi.complex = comp
            ddi_list.append(ddi)
        return ddi_list

    async def _run_docking(self):
        receptor = self._selected_receptor
        ligands = [item.complex for item in self._selected_ligands]
        site = Vector3(float(self._LocXInput.input_text), float(self._LocYInput.input_text), float(self._LocZInput.input_text))

        if not receptor or not ligands:
            if self._selected_site is None:
                Logs.warning("Trying to run docking without having one receptor, one site and at least one ligand selected")
                return

        self.show_loading(True)
        await self._plugin.run_docking(self._selected_receptor, ligands, site, self.get_params())
        self.show_loading(False)

    def make_plugin_usable(self, state=True):
        self._run_button.unusable = (not state) | self.refresh_run_btn_unusable(update=False)
        self._plugin.update_content(self._run_button)

    def show_loading(self, show=False):
        if show:
            self.ln_run_button.enabled = False
            self.ln_loading_bar.enabled = True
        else:
            self.ln_run_button.enabled = True
            self.ln_loading_bar.enabled = False
        self._plugin.update_menu(self._menu)

    def refresh_run_btn_unusable(self, update=True, after=False):
        site_requirement_met = self._selected_site is not None or not self._plugin._calculations.requires_site
        if self._selected_receptor is not None and len(self._selected_ligands) > 0 and site_requirement_met and not after:
            self._run_button.text.value.unusable = "Running..."
            self._run_button.text.size = 0.35
            self._run_button.unusable = False
        elif self._selected_receptor is not None and len(self._selected_ligands) > 0 and site_requirement_met and after:
            self._run_button.text.value.unusable = "Run"
            self._run_button.text.size = 0.35
            self._run_button.unusable = False
        else:
            self._run_button.text.value.unusable = "Please Select Complexes"
            self._run_button.text.size = 0.25
            self._run_button.unusable = True
        if update:
            self._plugin.update_content(self._run_button)

        return self._run_button.unusable

    def change_complex_list(self, complex_list):
        """When complex list is updated, repopulate dropdowns."""
        ligand_list = []
        receptor_list = []
        site_list = []

        ligand_list = self.create_complex_dropdown_items(complex_list)
        receptor_list = self.create_complex_dropdown_items(complex_list)
        site_list = self.create_complex_dropdown_items(complex_list)

        self.dd_ligands.items = ligand_list
        self.dd_receptor.items = receptor_list
        self.dd_site.items = site_list

        # Ligands should allow multiple selections
        for item in self.dd_ligands.items:
            item.close_on_selected = False

        # Reselect previously selected ligands
        for lig_item in self._selected_ligands:
            dd_item = next((item for item in self.dd_ligands.items if item.complex.index == lig_item.complex.index), None)
            if dd_item:
                dd_item.selected = True

        if not any([item.selected for item in self.dd_ligands.items]):
            self.dd_ligands.use_permanent_title = True
            self.dd_ligands.permanent_title = "None"
            self._selected_ligands = []

        # Reselect previously selected receptor
        if self._selected_receptor:
            new_receptor_item = next((item for item in self.dd_receptor.items if item.complex.index == self._selected_receptor.index), None)
            if new_receptor_item:
                new_receptor_item.selected = True

        if not any([item.selected for item in self.dd_receptor.items]):
            self.dd_receptor.use_permanent_title = True
            self.dd_receptor.permanent_title = "None"
            self._selected_receptor = None

        # Reselect previously selected site.
        if self._selected_site:
            new_site_item = next(
                (item for item in self.dd_site.items if item.complex.index == self._selected_site.complex.index), None)
            if new_site_item:
                new_site_item.selected = True

        if not any([item.selected for item in self.dd_site.items]):
            self.dd_site.use_permanent_title = True
            self.dd_site.permanent_title = "None"
            self._selected_site = None

        self.refresh_run_btn_unusable(update=False)
        self._plugin.update_menu(self._menu)

    def handle_ligand_selected(self, dropdown, item):
        self.multi_select_dropdown(dropdown, item)
        self._selected_ligands = dropdown._selected_items
        if not self._selected_ligands:
            self._ligand_txt._text_value = "Ligand"
        else:
            label = ', '.join([item.complex.full_name for item in self._selected_ligands])
            self._ligand_txt._text_value = label if len(label) <= 4 else label[:8] + '...'
        self.update_icons()
        self.refresh_run_btn_unusable(update=False)
        self._plugin.update_menu(self._menu)

    def handle_receptor_selected(self, dropdown, item):
        # If selected complex was previously selected, we are actually unselecting it.
        unselecting_complex = self._selected_receptor and item.complex.index == self._selected_receptor.index
        self._selected_receptor = None if unselecting_complex else item.complex

        if self._selected_receptor:
            self.dd_receptor.use_permanent_title = False
            self._receptor_txt._text_value = item.complex.full_name if len(item.complex.full_name) <= 4 else item.complex.full_name[:8] + '...'
        else:
            self._receptor_txt._text_value = "Receptor"
            self.dd_receptor.use_permanent_title = True
            self.dd_receptor.permanent_title = "None"
        self.update_icons()
        self.refresh_run_btn_unusable(update=False)
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
        self.make_plugin_usable()
        self._plugin.update_menu(self._menu)

    async def draw_site_sphere(self, comp, radius):
        Logs.debug('Drawing site sphere.')
        if hasattr(self, 'site_sphere'):
            Shape.destroy(self.site_sphere)
        self.site_sphere = Sphere()
        self.site_sphere.color = nanome.util.Color(100, 100, 100, 120)
        self.site_sphere.radius = radius
        anchor = self.site_sphere.anchors[0]
        anchor.anchor_type = nanome.util.enums.ShapeAnchorType.Complex
        anchor.target = comp.index
        complex_center = self.get_center(comp)
        anchor.local_offset = complex_center
        await Shape.upload(self.site_sphere)
        return self.site_sphere

    @async_callback
    async def handle_site_selected(self, dropdown, item):
        # If site was previously selected, we are actually unselecting it.
        unselecting_site = self._selected_site and item.complex.index == self._selected_site.complex.index
        self._selected_site = None if unselecting_site else item

        if self._selected_site:
            self.dd_site.use_permanent_title = False
            # Draw sphere indicating the site
            radius = self._slider.current_value
            comp = next(iter(await self._plugin.request_complexes([self._selected_site.complex.index])))
            asyncio.create_task(self.draw_site_sphere(comp, radius))
            complex_center = self.get_center(comp)
            self._LocXInput.input_text, self._LocYInput.input_text, self._LocZInput.input_text = [round(x, 2) for x in complex_center]
        else:
            self.dd_site.use_permanent_title = True
            self.dd_site.permanent_title = "None"
            item.selected = False
            self._LocXInput.input_text, self._LocYInput.input_text, self._LocZInput.input_text = '', '', ''
            if hasattr(self, 'site_sphere'):
                Shape.destroy(self.site_sphere)

        self.update_icons()
        self.refresh_run_btn_unusable(update=False)
        self._plugin.update_menu(self._menu)

    def update_icons(self):
        self._ligand_icon._file_path = ICONS['ligand_white' if self._selected_ligands else 'ligand_gray']
        self._receptor_icon._file_path = ICONS['receptor_white' if self._selected_receptor else 'receptor_gray']
        can_dock = self._selected_ligands and self._selected_receptor and self._selected_site
        self._check_arrow._file_path = ICONS['can_dock' if can_dock else 'cannot_dock']

    def build_menu(self):
        # panels
        root = self._menu.root
        self._docking_param_panel = root.find_node("LeftSide")
        self._score_panel = root.find_node("LeftSideScore")
        self._panel_separator = root.find_node("MiddleLine")

        # images
        self._check_arrow = root.find_node("CheckArrow").add_new_image(ICONS['cannot_dock'])
        self._ligand_icon = root.find_node("LigandIcon").add_new_image(ICONS['ligand_gray'])
        self._receptor_icon = root.find_node("ReceptorIcon").add_new_image(ICONS['receptor_gray'])

        slider_oval = root.find_node("SizeOval")
        slider_oval.add_new_image(file_path=ICONS['DarkOval'])

        refresh_icon = root.find_node("RefreshIcon")
        refresh_icon.add_new_image(file_path=ICONS['refresh'])

        # text
        self._txt2 = root.find_node("ModesInput").get_content()
        self._txt2.register_changed_callback(self.modes_changed)
        self._ligand_txt = root.find_node("LigandName").get_content()
        self._receptor_txt = root.find_node("ReceptorName").get_content()

        self._LocXInput = root.find_node("LocXInput").get_content()
        self._LocXInput.register_submitted_callback(partial(self.loc_submitted, 0))
        self._LocYInput = root.find_node("LocYInput").get_content()
        self._LocYInput.register_submitted_callback(partial(self.loc_submitted, 1))
        self._LocZInput = root.find_node("LocZInput").get_content()
        self._LocZInput.register_submitted_callback(partial(self.loc_submitted, 2))

        self.size_value_txt = root.find_node("SizeValue").get_content()

        align_btn = root.find_node("AlignButton").get_content()
        align_btn.register_pressed_callback(self.align_button_pressed_callback)
        align_btn.selected = True

        self._score_btn = root.find_node("ScoringButton").get_content()
        self._score_btn.register_pressed_callback(self.scoring_button_pressed_callback)

        close_score_btn = root.find_node("CloseScoreButton").get_content()
        close_score_btn.register_pressed_callback(self.close_score_pressed_callback)

        self._run_button.register_pressed_callback(self.run_button_pressed_callback)
        self._run_button.enabled = False
        self.refresh_run_btn_unusable(update=False)

        pose_sub_btn = root.find_node("PoseSub").get_content()
        pose_sub_btn.register_pressed_callback(self.pose_subbed_callback)
        pose_add_btn = root.find_node("PoseAdd").get_content()
        pose_add_btn.register_pressed_callback(self.pose_added_callback)

        location_refresh_btn = root.find_node("LocationRefresh").get_content()
        location_refresh_btn.register_pressed_callback(self.loc_refresh_pressed_callback)

        # dropdown
        self.dd_ligands.use_permanent_title = True
        self.dd_ligands.permanent_title = "None"
        self.dd_receptor.use_permanent_title = True
        self.dd_receptor.permanent_title = "None"
        self.dd_site.use_permanent_title = True
        self.dd_site.permanent_title = "None"

        # Ligands should allow multiple selections
        for item in self.dd_ligands.items:
            item.close_on_selected = False
        self.dd_ligands.register_item_clicked_callback(self.handle_ligand_selected)
        self.dd_receptor.register_item_clicked_callback(self.handle_receptor_selected)
        self.dd_site.register_item_clicked_callback(self.handle_site_selected)

        # slider
        self._slider = root.find_node("Slider").get_content()
        self._slider.register_released_callback(self.slider_released_callback)

        self._menu.enabled = True
        self._plugin.update_menu(self._menu)
        Logs.debug("Constructed plugin menu")

    @staticmethod
    def get_center(complex):
        """Calculate the center of a complex."""
        inf = float('inf')
        min_pos = Vector3(inf, inf, inf)
        max_pos = Vector3(-inf, -inf, -inf)

        for atom in complex.atoms:
            min_pos.x = min(min_pos.x, atom.position.x)
            min_pos.y = min(min_pos.y, atom.position.y)
            min_pos.z = min(min_pos.z, atom.position.z)
            max_pos.x = max(max_pos.x, atom.position.x)
            max_pos.y = max(max_pos.y, atom.position.y)
            max_pos.z = max(max_pos.z, atom.position.z)

        return (min_pos + max_pos) * 0.5

    @async_callback
    async def run_button_pressed_callback(self, button):
        if self._scoring:
            self._docking_param_panel.enabled = False
            self._score_panel.enabled = True
            self._score_list.items = []
            self._plugin.update_menu(self._menu)
        await self._run_docking()

    def modes_changed(self, input):
        try:
            self._modes = int(input.input_text)
        except:
            self._modes = 9

        if self._modes <= 0:
            self._modes = 1
            self._txt2.input_text = self._modes
            self._plugin.update_content(self._txt2)

    def align_button_pressed_callback(self, button):
        self._align = not self._align
        button.selected = self._align
        self._plugin.update_content(button)

    def close_score_pressed_callback(self, button):
        self._docking_param_panel.enabled = True
        self._score_panel.enabled = False
        self._plugin.update_menu(self._menu)

    def scoring_button_pressed_callback(self, button):
        self._scoring = not self._scoring
        button.selected = self._scoring
        self._plugin.update_content(button)

    @async_callback
    async def slider_released_callback(self, slider):
        slider.current_value = round(slider.current_value)
        self._plugin.update_content(slider)
        self._autobox = slider.current_value
        self.size_value_txt.text_value = str(slider.current_value)

        if hasattr(self, 'site_sphere'):
            # Update site sphere to show new radius
            self.site_sphere.radius = slider.current_value
            await Shape.upload(self.site_sphere)
        self._plugin.update_content(self.size_value_txt)

    def loc_submitted(index, self, text_input):
        try:
            float(text_input.input_text)
            self._selected_site.complex.position[index] = float(text_input.input_text)
        except:
            Logs.debug("Input is not a float")
        self._plugin.update_structures_shallow([self._selected_site.complex])

    def pose_added_callback(self, button):
        self._modes += 1
        self._txt2.input_text = self._modes
        self._plugin.update_content(self._txt2)

    def pose_subbed_callback(self, button):
        self._modes -= 1
        if self._modes <= 0:
            self._modes = 1
        self._txt2.input_text = self._modes
        self._plugin.update_content(self._txt2)

    def loc_refresh_pressed_callback(self, button):
        def update_site_loc(complexes_list):
            for complex in complexes_list:
                if complex.index == self._selected_site.complex.index:
                    self._selected_site.complex = complex
                    self._LocXInput.input_text, self._LocYInput.input_text, self._LocZInput.input_text = [round(x, 2) for x in complex.position]
                    self._plugin.update_menu(self._menu)

        if not self._selected_site:
            Logs.debug("No Site Selected")
            self._LocXInput.input_text, self._LocYInput.input_text, self._LocZInput.input_text = '', '', ''
            self._plugin.update_menu(self._menu)
        else:
            Logs.debug("Update the site location")
            self._plugin.request_complexes([self._selected_site.complex.index], update_site_loc)

    @property
    def dd_ligands(self):
        return self._menu.root.find_node("LigandDropdown").get_content()

    @property
    def dd_receptor(self):
        return self._menu.root.find_node("ReceptorDropdown").get_content()

    @property
    def dd_site(self):
        return self._menu.root.find_node("SiteDropdown").get_content()

    def multi_select_dropdown(self, dropdown, item):
        if not hasattr(dropdown, '_selected_items'):
            dropdown._selected_items = []

        selected_items = dropdown._selected_items
        if item not in selected_items:
            selected_items.append(item)
        else:
            selected_items.remove(item)
            item.selected = False

        for ddi in selected_items:
            ddi.selected = True

        dropdown.use_permanent_title = True
        permanent_title = ','.join([ddi.name for ddi in selected_items]) if selected_items else 'None'
        dropdown.permanent_title = permanent_title

        dropdown.use_permanent_title = len(selected_items) > 1
        self._plugin.update_content(dropdown)


class SettingsMenu:

    def __init__(self, plugin):
        self._plugin = plugin
        self._menu = nanome.ui.Menu.io.from_json(os.path.join(BASE_DIR, 'jsons', '_docking_setting_new.json'))
        self._exhaustiveness = 10

        menu_root = self._menu.root
        self.setting_slider_oval = menu_root.find_node("ExhaustOval")
        self.setting_slider_oval.add_new_image(file_path=ICONS['DarkOval'])
        self._exhaustiveness_txt = menu_root.find_node("ExhaustValue").get_content()
        self._display_score_btn = menu_root.find_node("VisualScoresButton").get_content()
        self._exhaust_slider = menu_root.find_node("ExhaustSlider").get_content()

        self._exhaust_slider.register_released_callback(self.exhaust_slider_released_callback)
        self._exhaust_slider.current_value = self._exhaustiveness
        self._display_score_btn.register_pressed_callback(self.visual_scores_button_pressed_callback)

    def exhaust_slider_released_callback(self, slider):
        slider.current_value = round(slider.current_value)
        self._plugin.update_content(slider)
        self._exhaustiveness = slider.current_value
        self._exhaustiveness_txt.text_value = str(self._exhaustiveness)
        self._plugin.update_content(self._exhaustiveness_txt)

    def visual_scores_button_pressed_callback(self, button):
        self._visual_scores = not self._visual_scores
        button.selected = self._visual_scores
        self._plugin.update_content(button)
