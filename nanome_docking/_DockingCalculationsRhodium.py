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
        self.__docking_running = False

    def start_docking(self, receptor, ligands, site, params):
        self.initialize()
        self._ligands = ligands
        self._receptor = receptor
        self._params = params
        self.prepare_receptor(receptor)

    def update(self):
        if self.__docking_running == False:
            return

        self.__process.communicate()
        if self.__process.poll() != None:
            self.__docking_running = False
            self._docking_finished()

    def initialize(self):
        self._protein_input = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb")
        self._protein_converted_input = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb")
        self._ligands_input = tempfile.NamedTemporaryFile(delete=False, suffix=".sdf")
        self._ligands_converted_input = tempfile.NamedTemporaryFile(delete=False, suffix=".sdf")
        self._docking_output = tempfile.NamedTemporaryFile()
        self._docking_output.close()

    def clean_files(self, docking_result):
        self._protein_input.close()
        self._protein_converted_input.close()
        self._ligands_input.close()
        self._ligands_converted_input.close()
        os.remove(self._protein_input.name)
        os.remove(self._protein_converted_input.name)
        os.remove(self._ligands_input.name)
        os.remove(self._ligands_converted_input.name)
        os.remove(self._docking_output.name + ".csv")
        for entry in docking_result:
            os.chmod(entry[0], 0o777)
            os.remove(entry[0])
        pass

    def prepare_receptor(self, receptor):
        residues_to_remove = []
        for residue in receptor.residues:
            if residue.name == "HOH":
                residues_to_remove.append(residue)

        for residue in residues_to_remove:
            residue.parent.remove_residue(residue)

        receptor.io.to_pdb(self._protein_input.name)

        proc = Process()
        proc.executable_path = 'obabel'
        proc.args = ['-ipdb', self._protein_input.name,
            '-opdb', '-O' + self._protein_converted_input.name, '-h']
        proc.on_done = self.receptor_ready
        proc.start()

    def receptor_ready(self, return_value):
        self._ligands.io.to_sdf(self._ligands_input.name)
        self._start_conversion()

    def _start_conversion(self):
        proc = Process()
        proc.executable_path = 'obabel'
        proc.args = ['-isdf', self._ligands_input.name,
            '-osdf', '-O' + self._ligands_converted_input.name]
        proc.on_done = self._conversion_finished
        proc.start()

    def _conversion_finished(self, return_value):
        self._start_docking()

    def _start_docking(self):
        args = [RHODIUM_PATH, self._protein_converted_input.name, self._ligands_converted_input.name,
            '--outfile', self._docking_output.name,
            '--refine', str(self._params['poses']),
            '--resolution', str(self._params['grid_resolution']),
            '--refine', str(self._params['poses']),
            '--nr', str(self._params['rotamers'])]

        if self._params['ignore_hetatoms']:
            args += ['--ignore_pdb_hetatm']

        Logs.debug("Start Rhodium:", args)

        self._start_timer = timer()
        try:
            self.__process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=os.path.dirname(__file__))
        except:
            nanome.util.Logs.error("Couldn't execute Rhodium, please check if executable is in the plugin folder and has permissions. Path:", exe_path, traceback.format_exc())
            self._plugin.make_plugin_usable()
            self._plugin.send_notification(NotificationTypes.error, "Docking error, check plugin")
            return

        self._plugin.send_notification(NotificationTypes.message, "Docking started")
        self.__docking_running = True

    def _on_docking_output(self, output):
        Logs.debug("Docking output:", output)
        pass

    def _on_docking_error(self, error):
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
            complex = nanome.structure.Complex.io.from_pdb(path=entry[0])
            molecule = next(complex.molecules)
            molecule._associated['score'] = entry[1]
            molecule._associated['affinity'] = entry[2]
            molecule._associated['pose_population'] = entry[3]
            docked_ligands.add_molecule(molecule)

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


    def _docking_finished(self):
        end = timer()
        nanome.util.Logs.debug("Docking Finished in", end - self._start_timer, "seconds")

        docking_output = self.read_csv()
        self.assemble_result(docking_output)

        self.clean_files(docking_output)
        self._plugin.make_plugin_usable()
        self._plugin.send_notification(NotificationTypes.success, "Docking finished")
