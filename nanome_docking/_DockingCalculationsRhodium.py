import nanome
from nanome.util import Logs, Process
from nanome.util.enums import NotificationTypes

import os
import tempfile
import traceback
import subprocess
from timeit import default_timer as timer

RHODIUM_PATH = os.path.join(os.path.dirname(__file__), 'Rh_x64.exe')

class DockingCalculations():
    def __init__(self, plugin):
        self._plugin = plugin

    def start_docking(self, receptor, ligands, site, params):
        self.initialize()

        receptor.io.to_pdb(self._protein_input.name)
        nanome.util.Logs.debug("Saved PDB", self._protein_input.name)
        ligands.io.to_sdf(self._ligands_input.name)
        nanome.util.Logs.debug("Saved SDF", self._ligands_input.name)

        self._receptor = receptor
        self._params = params

        # Start docking process
        self._start_docking()

    def initialize(self):
        self._protein_input = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb")
        self._ligands_input = tempfile.NamedTemporaryFile(delete=False, suffix=".sdf")
        self._docking_output = tempfile.NamedTemporaryFile()
        self._docking_output.close()

    def clean_files(self, docking_result):
        os.remove(self._protein_input.name)
        os.remove(self._ligands_input.name)
        os.remove(self._docking_output.name + ".csv")
        for entry in docking_result:
            os.remove(entry[0])
        pass

    def update(self):
        return

    def _start_docking(self):
        self.__proc = Process()
        self.__proc.executable_path = RHODIUM_PATH
        self.__proc.args = [self._protein_input.name, self._ligands_input.name,
            '--outfile', self._docking_output.name,
            '--refine', str(self._params['poses']),
            '--resolution', str(self._params['grid_resolution']),
            '--refine', str(self._params['poses']),
            '--nr', str(self._params['rotamers'])]
        self.__proc.output_text = True

        if self._params['ignore_hetatoms']:
            args += ['--ignore_pdb_hetatm']

        nanome.util.Logs.debug("Run Rhodium")
        self._start_timer = timer()
        try:
            self.__proc.on_error = self._docking_error
            self.__proc.on_done = self._docking_finished
            self.__proc.start()
        except:
            nanome.util.Logs.error("Couldn't execute Rhodium, please check if executable is in the plugin folder and has permissions. Path:", exe_path, traceback.format_exc())
            self._plugin.make_plugin_usable()
            self._plugin.send_notification(NotificationTypes.error, "Docking error, check plugin")
            return

        self._plugin.send_notification(NotificationTypes.message, "Docking started")

    def _docking_error(self, error):
        Logs.error("Docking error:", error)

    def read_csv(self):
        result = []
        with open(self._docking_output.name + ".csv", 'r') as file:
            lines = file.readlines()
            for line in lines:
                if line.startswith('seed') == False:
                    continue
                line_split = line.split(',')
                result.append([line_split[5].strip(),
                    line_split[13].strip(),
                    line_split[41].strip(),
                    line_split[43].strip()])
        return result

    def assemble_result(self, docking_output):
        docked_ligands = nanome.api.structure.Complex()
        for entry in docking_output:
            complex = nanome.structure.Complex.io.from_sdf(path=entry[0])
            docked_ligands.add_molecule(complex.molecules[0])

        def bonds_added(complex_arr):
            docked_ligands_bonded = complex_arr[0]
            docked_ligands_bonded.molecular.name = "Docking"
            docked_ligands_bonded.rendering.visible = True
            if self._params['align'] == True:
                docked_ligands_bonded.transform.position = self._receptor.transform.position
                docked_ligands_bonded.transform.rotation = self._receptor.transform.rotation
                docked_ligands_bonded.rendering.boxed = True

            nanome.util.Logs.debug("Update workspace")
            self._plugin.add_result_to_workspace([docked_ligands_bonded])

        self._plugin.add_bonds([docked_ligands], bonds_added)


    def _docking_finished(self, return_value):
        end = timer()
        nanome.util.Logs.debug("Docking Finished in", end - self._start_timer, "seconds")

        docking_output = self.read_csv()
        self.assemble_result(docking_output)

        self.clean_files(docking_output)
        self._plugin.make_plugin_usable()
        self._plugin.send_notification(NotificationTypes.success, "Docking finished")
