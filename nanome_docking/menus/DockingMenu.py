import os
import nanome
from nanome.util import Logs, async_callback
from nanome.api.ui import DropdownItem
from nanome.api.shapes import Sphere, Shape
from nanome_docking.utils import get_complex_center


BASE_DIR = os.path.dirname(__file__)
ICONS_DIR = os.path.join(BASE_DIR, 'icons')
ICONS = {icon.rsplit('.')[0]: os.path.join(ICONS_DIR, icon) for icon in os.listdir(ICONS_DIR)}

SETTINGS_JSON_PATH = os.path.join(BASE_DIR, 'jsons', 'docking_settings.json')
MENU_JSON_PATH = os.path.join(BASE_DIR, 'jsons', 'docking_menu.json')


class DockingMenu():

    def __init__(self, docking_plugin, algorithm=None):
        algorithm = algorithm or ''
        self._plugin = docking_plugin

        self._menu = nanome.ui.Menu.io.from_json(MENU_JSON_PATH)
        self._menu.title = f'{algorithm} Docking'
    
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


        # Run button
        self.ln_run_button = self._menu.root.find_node("RunButton")
        self._run_button = self.ln_run_button.get_content()
        # loading bar
        self.ln_loading_bar = self._menu.root.find_node("LoadingBar")
        self.loading_bar = self.ln_loading_bar.get_content()


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
        ligands = self._selected_ligands
        site = None
        if self._selected_site:
            site = self._selected_site

        if not receptor or not ligands:
            if self._selected_site is None:
                Logs.warning("Trying to run docking without having one receptor, one site and at least one ligand selected")
                return

        self.loading_bar.percentage = 0
        self.enable_loading_bar()
        self.make_plugin_usable(False)
        await self._plugin.run_docking(self._selected_receptor, ligands, site, self.get_params())
        self.enable_loading_bar(False)

    def make_plugin_usable(self, state=True):
        self._run_button.unusable = (not state) or self.refresh_run_btn_unusable(update=False)
        self._plugin.update_content(self._run_button)

    def refresh_run_btn_unusable(self, update=True, after=False):
        site_requirement_met = self._selected_site is not None
        if self._selected_receptor is not None and len(self._selected_ligands) > 0 and site_requirement_met and not after:
            self._run_button.text.value.unusable = "Running..."
            self._run_button.text.size = 0.35
            self._run_button.unusable = False
        elif self._selected_receptor is not None and len(self._selected_ligands) > 0 and site_requirement_met and after:
            self._run_button.text.value.unusable = "Run"
            self._run_button.text.size = 0.35
            self._run_button.unusable = False
        else:
            self._run_button.text.value.unusable = "Run"
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
        for comp in self._selected_ligands:
            dd_item = next((item for item in self.dd_ligands.items if item.complex.index == comp.index), None)
            if dd_item:
                dd_item.selected = True

        if not any([item.selected for item in self.dd_ligands.items]):
            self.dd_ligands.use_permanent_title = True
            self.dd_ligands.permanent_title = "None"
            self._selected_ligands = []

        # Reselect previously selected receptor
        if self._selected_receptor:
            new_receptor_item = next(
                (item for item in self.dd_receptor.items if item.complex.index == self._selected_receptor.index), None)
            if new_receptor_item:
                new_receptor_item.selected = True

        if not any([item.selected for item in self.dd_receptor.items]):
            self.dd_receptor.use_permanent_title = True
            self.dd_receptor.permanent_title = "None"
            self._selected_receptor = None

        # Reselect previously selected site.
        if self._selected_site:
            new_site_item = next(
                (item for item in self.dd_site.items if item.complex.index == self._selected_site.index), None)
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

        dropdown.use_permanent_title = len(dropdown._selected_items) > 1
        if not dropdown._selected_items:
            ligand_text = "Ligands"
            dropdown.permanent_title = "None"
        else:
            complex_names = ','.join([ddi.complex.full_name for ddi in dropdown._selected_items])
            ligand_text = complex_names if len(complex_names) <= 8 else f'{complex_names[:7]}...'
            dropdown.permanent_title = complex_names

        self._ligand_txt._text_value = ligand_text
        self._selected_ligands = [ddi.complex for ddi in dropdown._selected_items]

        self.update_icons()
        self.refresh_run_btn_unusable(update=False)
        self._plugin.update_menu(self._menu)

    def handle_receptor_selected(self, dropdown, item):
        # If selected complex was previously selected, we are actually unselecting it.
        unselecting_complex = self._selected_receptor and item.complex.index == self._selected_receptor.index
        self._selected_receptor = None if unselecting_complex else item.complex

        if self._selected_receptor:
            dropdown.use_permanent_title = False
            self._receptor_txt._text_value = item.complex.full_name if len(item.complex.full_name) <= 4 else item.complex.full_name[:8] + '...'
        else:
            self._receptor_txt._text_value = "Receptor"
            dropdown.use_permanent_title = True
            dropdown.permanent_title = "None"
        self.update_icons()
        self.refresh_run_btn_unusable(update=False)
        self._plugin.update_menu(self._menu)

    @async_callback
    async def draw_site_sphere(self, site_complex, radius):
        """Draw sphere at origin with provided radius.

        :arg origin: Vector3, center of the site sphere
        :arg radius: int, radius of sphere.
        """
        Logs.debug('Drawing site sphere.')
        if not hasattr(self, 'site_sphere'):
            self.site_sphere = Sphere()
            self.site_sphere.color = nanome.util.Color(0, 100, 0, 120)
        self.site_sphere.radius = radius
        anchor = self.site_sphere.anchors[0]
        anchor.anchor_type = nanome.util.enums.ShapeAnchorType.Workspace
        site_center = get_complex_center(site_complex)
        anchor.local_offset = site_complex.get_complex_to_workspace_matrix() * site_center
        await Shape.upload(self.site_sphere)
        return self.site_sphere

    @async_callback
    async def handle_site_selected(self, dropdown, item):
        # If site was previously selected, we are actually unselecting it.
        unselecting_site = self._selected_site and item.complex.index == self._selected_site.index
        self._selected_site = None if unselecting_site else item.complex

        if self._selected_site:
            dropdown.use_permanent_title = False
            comp = next(iter(await self._plugin.request_complexes([self._selected_site.index])))
            # Draw sphere indicating the site
            radius = self._slider.current_value
            self.draw_site_sphere(comp, radius)
            complex_center = get_complex_center(comp)
            self._site_x.input_text, self._site_y.input_text, self._site_z.input_text = [
                round(x, 2) for x in complex_center]
        else:
            dropdown.use_permanent_title = True
            dropdown.permanent_title = "None"
            item.selected = False
            self._site_x.input_text, self._site_y.input_text, self._site_z.input_text = '', '', ''
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

    def build_menu(self, current_algorithm):
        self._menu.title = f'{current_algorithm.title()} Docking'
        # panels
        root = self._menu.root
        self._docking_param_panel = root.find_node("LeftSide")
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

        self._site_x = root.find_node("LocXInput").get_content()
        self._site_y = root.find_node("LocYInput").get_content()
        self._site_z = root.find_node("LocZInput").get_content()
        self._site_x.input_text, self._site_y.input_text, self._site_z.input_text = '', '', ''

        self._site_x.input_text
        self.size_value_txt = root.find_node("SizeValue").get_content()

        self._run_button.register_pressed_callback(self.run_button_pressed_callback)
        self._run_button.enabled = False
        self.refresh_run_btn_unusable(update=False)

        pose_sub_btn = root.find_node("PoseSub").get_content()
        pose_sub_btn.register_pressed_callback(self.pose_subbed_callback)
        pose_add_btn = root.find_node("PoseAdd").get_content()
        pose_add_btn.register_pressed_callback(self.pose_added_callback)

        location_refresh_btn = root.find_node("LocationRefresh").get_content()
        location_refresh_btn.register_pressed_callback(self.loc_refresh_pressed_callback)

        # dropdowns
        self.dd_ligands = self._menu.root.find_node("LigandDropdown").get_content()
        self.dd_receptor = self._menu.root.find_node("ReceptorDropdown").get_content()
        self.dd_site = self._menu.root.find_node("SiteDropdown").get_content()
        for dd in [self.dd_ligands, self.dd_receptor, self.dd_site]:
            dd.use_permanent_title = True
            dd.permanent_title = "None"

        # Ligands should allow multiple selections
        for item in self.dd_ligands.items:
            item.close_on_selected = False
        self.dd_ligands.register_item_clicked_callback(self.handle_ligand_selected)
        self.dd_receptor.register_item_clicked_callback(self.handle_receptor_selected)
        self.dd_site.register_item_clicked_callback(self.handle_site_selected)

        # slider
        self._slider = root.find_node("Slider").get_content()
        self._slider.register_released_callback(self.slider_released_callback)

        self._plugin.update_menu(self._menu)
        Logs.debug("Constructed plugin menu")

    def enable(self):
        self._menu.enabled = True
        self._plugin.update_menu(self._menu)

    @async_callback
    async def run_button_pressed_callback(self, button):
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

    @async_callback
    async def loc_refresh_pressed_callback(self, button):

        if not self._selected_site:
            Logs.debug("No Site Selected")
            self._site_x.input_text, self._site_y.input_text, self._site_z.input_text = '', '', ''
            self._plugin.update_menu(self._menu)
        else:
            Logs.debug("Update the site location")
            comp = self._selected_site
            comp = (await self._plugin.request_complexes([comp.index]))[0]
            self._site_x.input_text, self._site_y.input_text, self._site_z.input_text = [
                round(x, 2) for x in comp.position]
            radius = self._slider.current_value
            self.draw_site_sphere(comp, radius)
            self._plugin.update_menu(self._menu)

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

    def enable_loading_bar(self, enabled=True):
        self.ln_loading_bar.enabled = enabled
        self._plugin.update_node(self.ln_loading_bar)

    def update_loading_bar(self, current, total):
        self.loading_bar.percentage = current / total
        self._plugin.update_content(self.loading_bar)


class SettingsMenu:

    def __init__(self, plugin):
        self.plugin = plugin
        self._menu = nanome.ui.Menu.io.from_json(SETTINGS_JSON_PATH)
        self._menu.index = 1
        self._exhaustiveness = 10

        menu_root = self._menu.root
        self.setting_slider_oval = menu_root.find_node("ExhaustOval")
        self.setting_slider_oval.add_new_image(file_path=ICONS['DarkOval'])
        self._exhaustiveness_txt = menu_root.find_node("ExhaustValue").get_content()
        self._exhaust_slider = menu_root.find_node("ExhaustSlider").get_content()
        self._visual_scores = False

        self.dd_algorithm = menu_root.find_node("dd_algorithm").get_content()
        self.dd_algorithm.register_item_clicked_callback(self.on_algorithm_selected)

        self._exhaust_slider.register_released_callback(self.exhaust_slider_released_callback)
        self._exhaust_slider.current_value = self._exhaustiveness

    def enable(self):
        self._menu.enabled = True
        self.plugin.update_menu(self._menu)

    def exhaust_slider_released_callback(self, slider):
        slider.current_value = round(slider.current_value)
        self.plugin.update_content(slider)
        self._exhaustiveness = slider.current_value
        self._exhaustiveness_txt.text_value = str(self._exhaustiveness)
        self.plugin.update_content(self._exhaustiveness_txt)

    @property
    def current_algorithm(self):
        selected_dd_item = next(item for item in self.dd_algorithm.items if item.selected)
        algorithm = selected_dd_item.name.lower()
        return algorithm

    def on_algorithm_selected(self, dropdown, item):
        self.plugin.change_algorithm(item.name)
        

        

