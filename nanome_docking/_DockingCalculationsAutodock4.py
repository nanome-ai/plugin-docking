import nanome
import os
import shutil
import subprocess
import tempfile
import re
from timeit import default_timer as timer

# TMP
from nanome._internal._structure._io._pdb.save import Options as PDBOptions
from nanome._internal._structure._io._sdf.save import Options as SDFOptions

from .ComplexUtils import ComplexUtils

class DockingCalculations():
    def __init__(self, plugin):
        self._plugin = plugin
        self._preparation_pending = False
        self._grid_pending = False
        self._docking_pending = False
        self._bond_pending = False
        self._running = False
        self._pdb_options = PDBOptions()
        self._pdb_options.write_bonds = True
        self._sdf_options = SDFOptions()
        self._sdf_options.write_bonds = True
        self.requires_site = False

    def initialize(self):
        # TODO: Read and write in a folder unique per plugin instance
        self.temp_dir = tempfile.TemporaryDirectory()
        self._protein_input = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb", dir=self.temp_dir.name)
        self._protein_input_converted = tempfile.NamedTemporaryFile(delete=False, suffix=".pdbqt", dir=self.temp_dir.name)
        self._ligands_input = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb", dir=self.temp_dir.name)
        self._ligands_input_converted = tempfile.NamedTemporaryFile(delete=False, suffix=".pdbqt", dir=self.temp_dir.name)
        self._ligands_output = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb", dir=self.temp_dir.name)
        self._bond_output = tempfile.NamedTemporaryFile(delete=False, suffix=".sdf", dir=self.temp_dir.name)
        self._autogrid_input = tempfile.NamedTemporaryFile(delete=False, suffix=".gpf", dir=self.temp_dir.name)
        self._autodock_input = tempfile.NamedTemporaryFile(delete=False, suffix=".dpf", dir=self.temp_dir.name)
        self._autogrid_log = tempfile.NamedTemporaryFile(delete=False, suffix=".glg", dir=self.temp_dir.name)
        self._autodock_log = tempfile.NamedTemporaryFile(delete=False, suffix=".dlf", dir=self.temp_dir.name)

    def start_docking(self, receptor, ligands, site, exhaustiveness, modes, align, replace, scoring, visual_scores, autobox):
        self.initialize()
        # Save all input files
        receptor.io.to_pdb(self._protein_input.name, self._pdb_options)
        nanome.util.Logs.debug("Saved PDB", self._protein_input.name)
        ligands.io.to_pdb(self._ligands_input.name, self._pdb_options)
        nanome.util.Logs.debug("Saved PDB", self._ligands_input.name)

        self._receptor = receptor
        self._combined_ligands = ComplexUtils.combine_ligands(receptor, ligands)
        self._site = site
        self._align = align
        self._replace = replace
        self._visual_scores = visual_scores

        # Start docking process
        self._running = False
        self._preparation_pending = True
        self._parameters_preparation_pending = True
        self._grid_pending = True
        self._docking_pending = True
        self._bond_pending = True

    def update(self):
        if self._bond_pending == False:
            return

        if self._running == False:
            if self._preparation_pending == True:
                self._start_preparation()
            elif self._parameters_preparation_pending == True:
                self._start_parameters_preparation()
            elif self._grid_pending == True:
                self._start_grid()
            elif self._docking_pending == True:
                self._start_docking()
            elif self._bond_pending == True:
                self._start_bonds()
        elif self._preparation_pending:
            if self._check_preparation():
                self._preparation_finished()
        elif self._parameters_preparation_pending:
            if self._check_parameters_preparation():
                self._parameters_preparation_finished()
        elif self._grid_pending:
            if self._check_grid():
                self._grid_finished()
        elif self._docking_pending:
            if self._check_docking():
                self._docking_finished()
        elif self._bond_pending:
            if self._check_bonds():
                self._bonds_finished()

    def _check_process_error(self, process, check_only_errors = False):
        (results, errors) = process.communicate()

        try:
            if len(errors) == 0:
                if check_only_errors == True:
                    return True
                for line in results.splitlines():
                    nanome.util.Logs.debug(line.decode("utf-8"))
            else:
                error = False
                for line in errors.splitlines():
                    str = line.decode("utf-8")
                    if "warning" in str.lower():
                        nanome.util.Logs.warning(str)
                    elif str.trim() != "":
                        nanome.util.Logs.error(str)
                        error = True
                if error == True:
                    self._preparation_pending = False
                    self._parameters_preparation_pending = False
                    self._grid_pending = False
                    self._docking_pending = False
                    self._bond_pending = False
                    self._running = False
                    self._plugin.make_plugin_usable()
                    nanome.util.Logs.error("Error. Abort docking")
                    return True
        except:
            pass
        return False

    # Preparation of lig and receptor files

    def _start_preparation(self):
        # Awful situation here
        lig_args = ['python', os.path.join(os.path.dirname(__file__), 'prepare_ligand4.py'), '-l', self._ligands_input.name, '-o', self._ligands_input_converted.name]
        rec_args = ['python', os.path.join(os.path.dirname(__file__), 'prepare_receptor4.py'), '-r', self._protein_input.name, '-o', self._protein_input_converted.name]

        nanome.util.Logs.debug("Prepare ligand and receptor")
        self._start_timer = timer()
        self._lig_process = subprocess.Popen(lig_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=self.temp_dir.name)
        self._rec_process = subprocess.Popen(rec_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=self.temp_dir.name)
        self._running = True

    def _check_preparation(self):
        return self._lig_process.poll() != None or self._rec_process.poll() != None

    def _preparation_finished(self):
        end = timer()
        nanome.util.Logs.debug("Prepared ligand and receptor in", end - self._start_timer, "seconds")

        if self._check_process_error(self._lig_process) or self._check_process_error(self._rec_process):
            return

        self._running = False
        self._preparation_pending = False
        self._parameters_preparation_pending = True

    # Preparation of gpf and dpf files

    def _start_parameters_preparation(self):
        # Awful situation here
        grid_args = ['python', os.path.join(os.path.dirname(__file__), 'prepare_gpf4.py'), '-l', self._ligands_input_converted.name, '-r', self._protein_input_converted.name, '-o', self._autogrid_input.name]
        dock_args = ['python', os.path.join(os.path.dirname(__file__), 'prepare_dpf42.py'), '-l', self._ligands_input_converted.name, '-r', self._protein_input_converted.name, '-o', self._autodock_input.name]

        nanome.util.Logs.debug("Prepare grid and docking parameter files")
        self._start_timer = timer()
        self._grid_process = subprocess.Popen(grid_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=self.temp_dir.name)
        self._dock_process = subprocess.Popen(dock_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=self.temp_dir.name)
        self._running = True

    def _check_parameters_preparation(self):
        return self._grid_process.poll() != None or self._dock_process.poll() != None

    def _parameters_preparation_finished(self):
        end = timer()
        nanome.util.Logs.debug("Prepared grid and docking parameter files in", end - self._start_timer, "seconds")

        if self._check_process_error(self._grid_process) or self._check_process_error(self._dock_process):
            return

        f = open(self._autodock_input.name, 'a')
        f.write("write_all\n")
        f.close()

        self._running = False
        self._parameters_preparation_pending = False
        self._grid_pending = True

    # Autogrid

    def _start_grid(self):
        full_name = self._autogrid_input.name
        full_name_log = self._autogrid_log.name
        path = os.path.dirname(full_name)

        args = ['autogrid4', '-p', full_name, '-l', full_name_log]

        nanome.util.Logs.debug("Start Autogrid")
        self._start_timer = timer()
        self._autogrid_process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=path)
        self._running = True

    def _check_grid(self):
        return self._autogrid_process.poll() != None

    def _grid_finished(self):
        end = timer()
        nanome.util.Logs.debug("Ran Autogrid in", end - self._start_timer, "seconds. Logs: ", self._autogrid_log.name)

        if self._check_process_error(self._autogrid_process):
            return

        self._running = False
        self._grid_pending = False
        self._docking_pending = True

    # Autodock

    def _start_docking(self):
        full_name_input = self._autodock_input.name
        path = os.path.dirname(full_name_input)
        full_name_log = self._autodock_log.name

        args = ['autodock4', '-p', full_name_input, '-l', full_name_log]

        nanome.util.Logs.debug("Start Autodock")
        self._start_timer = timer()
        self._autodock_process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=path)
        self._running = True

    def _check_docking(self):
        return self._autodock_process.poll() != None

    def _docking_finished(self):
        end = timer()
        nanome.util.Logs.debug("Ran Autodock in", end - self._start_timer, "seconds. Logs:", self._autodock_log.name)

        if self._check_process_error(self._autodock_process):
            return

        # Conversion

        with open(self._autodock_log.name) as origin_file, open(self._ligands_output.name, 'w') as destination_file:
            print("ligand output file:", self._ligands_output.name)

            for line in origin_file:
                if line.startswith('DOCKED'):
                    str = line[8:]
                    if str.startswith('ROOT') or str.startswith('ENDROOT'):
                        continue
                    str = str.replace('+', ' ')
                    if not str.endswith('\n'):
                        str += '\n'
                    str = re.sub("([0-9]{3}) HD *\n", r"\1 H \n", str)
                    str = re.sub("([0-9]{3}) N.*A *\n", r"\1 N \n", str)
                    str = re.sub("([0-9]{3}) O.*A *\n", r"\1 O \n", str)
                    str = re.sub("([0-9]{3}) S.*A *\n", r"\1 S \n", str)
                    str = re.sub("([0-9]{3}) A.* *\n", r"\1 C \n", str)
                    str = re.sub("([0-9]{3} [A-G][a-g]) \n", r"\1\n", str)
                    str = re.sub("([0-9]{3}) ([A-G][a-g])\n", r"\1\2\n", str)
                    destination_file.write(str)

        self._running = False
        self._docking_pending = False
        self._bond_pending = True

    # Add Bonds

    def _start_bonds(self):
        nanome.util.Logs.debug("Start Bonds")
        self._start_timer = timer()
        args = ['obabel', '-ipdb', self._ligands_output.name, '-osdf', '-O' + self._bond_output.name]
        self._obabel_process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self._running = True

    def _check_bonds(self):
        poll = self._obabel_process.poll()
        if poll == None:
            self._obabel_process.communicate()
        return poll != None

    def _bonds_finished(self):
        end = timer()
        nanome.util.Logs.debug("Ran OpenBabel in", end - self._start_timer, "seconds")

        # if self._check_process_error(self._obabel_process, check_only_errors=True):
        #     return

        self._running = False
        self._bond_pending = False

        docked_ligands = nanome.structure.Complex.io.from_sdf(path=self._bond_output.name)
        nanome.util.Logs.debug("Read SDF", self._bond_output.name)

        docked_ligands.name = self._combined_ligands.full_name + " (Docked)"
        docked_ligands.visible = True
        if self._align == True:
            docked_ligands.transform.position = self._receptor.transform.position
            docked_ligands.transform.rotation = self._receptor.transform.rotation

        nanome.util.Logs.debug("Update workspace")
        # TODO: verify this shouldn't be here anymore (test)
        # self._plugin.make_plugin_usable()
        self._plugin.add_result_to_workspace([docked_ligands])

        # shutil.rmtree(self.temp_dir.name)