import nanome
import gzip
import itertools
import operator
import os
import subprocess
import tempfile
import re
from timeit import default_timer as timer

# TMP
from nanome._internal._structure._io._pdb.save import Options as PDBOptions
from nanome._internal._structure._io._sdf.save import Options as SDFOptions

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

    def initialize(self):
        # TODO: Read and write in a folder unique per plugin instance
        self._protein_input = os.getcwd() + ("/protein.pdb")
        self._protein_input_converted = os.getcwd() + ("/protein.pdbqt")
        self._ligand_input = os.getcwd() + ("/ligand.pdb")
        self._ligand_input_converted = os.getcwd() + ("/ligand.pdbqt")
        self._ligand_output = os.getcwd() + ("/output.pdb")
        self._bond_output = os.getcwd() + ("/output.sdf")
        self._autogrid_input = os.getcwd() + ("/grid.gpf")
        self._autodock_input = os.getcwd() + ("/dock.dpf")
        self._autodock_input = os.getcwd() + ("/dock.dpf")
        self._autodock_log = os.getcwd() + ("/dock_log.dlg")

    def start_docking(self, receptor, ligand, site, exhaustiveness, modes, align, replace):
        # Save all input files
        receptor.io.to_pdb(self._protein_input, self._pdb_options)
        nanome.util.Logs.debug("Saved PDB", self._protein_input)
        ligand.io.to_pdb(self._ligand_input, self._pdb_options)
        nanome.util.Logs.debug("Saved PDB", self._ligand_input)

        self._receptor = receptor
        self._site = site
        self._align = align
        self._replace = replace

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
                for line in errors.splitlines():
                    nanome.util.Logs.error(line.decode("utf-8"))
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
        lig_args = ['py', '-2.5', 'prepare_ligand4.py', '-l', self._ligand_input, '-o', self._ligand_input_converted]
        rec_args = ['py', '-2.5', 'prepare_receptor4.py', '-r', self._protein_input, '-o', self._protein_input_converted]
        
        nanome.util.Logs.debug("Prepare ligand and receptor")
        self._start_timer = timer()
        self._lig_process = subprocess.Popen(lig_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self._rec_process = subprocess.Popen(rec_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
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
        grid_args = ['py', '-2.5', 'prepare_gpf4.py', '-l', self._ligand_input_converted, '-r', self._protein_input_converted, '-o', self._autogrid_input]
        dock_args = ['py', '-2.5', 'prepare_dpf42.py', '-l', self._ligand_input_converted, '-r', self._protein_input_converted, '-o', self._autodock_input]
        
        nanome.util.Logs.debug("Prepare grid and docking parameter files")
        self._start_timer = timer()
        self._grid_process = subprocess.Popen(grid_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self._dock_process = subprocess.Popen(dock_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self._running = True

    def _check_parameters_preparation(self):
        return self._grid_process.poll() != None or self._dock_process.poll() != None

    def _parameters_preparation_finished(self):
        end = timer()
        nanome.util.Logs.debug("Prepared grid and docking parameter files in", end - self._start_timer, "seconds")

        if self._check_process_error(self._grid_process) or self._check_process_error(self._dock_process):
            return

        f = open(self._autodock_input, 'a')
        f.write("write_all\n")
        f.close()

        self._running = False
        self._parameters_preparation_pending = False
        self._grid_pending = True

    # Autogrid

    def _start_grid(self):
        args = ['autogrid4', '-p', self._autogrid_input, '-l', 'grid_log.glg']
        
        nanome.util.Logs.debug("Start Autogrid")
        self._start_timer = timer()
        self._autogrid_process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self._running = True

    def _check_grid(self):
        return self._autogrid_process.poll() != None

    def _grid_finished(self):
        end = timer()
        nanome.util.Logs.debug("Ran Autogrid in", end - self._start_timer, "seconds")

        if self._check_process_error(self._autogrid_process):
            return

        self._running = False
        self._grid_pending = False
        self._docking_pending = True

    # Autodock

    def _start_docking(self):
        args = ['autodock4', '-p', self._autodock_input, '-l', self._autodock_log]
        
        nanome.util.Logs.debug("Start Autodock")
        self._start_timer = timer()
        self._autodock_process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self._running = True

    def _check_docking(self):
        return self._autodock_process.poll() != None

    def _docking_finished(self):
        end = timer()
        nanome.util.Logs.debug("Ran Autodock in", end - self._start_timer, "seconds")

        if self._check_process_error(self._autodock_process):
            return

        # Conversion

        with open(self._autodock_log) as origin_file, open(self._ligand_output, 'w') as destination_file:
            for line in origin_file:
                if line.startswith('DOCKED'):
                    str = line[8:]
                    if str.startswith('ROOT') or str.startswith('ENDROOT'):
                        continue
                    idx = str.find('+')
                    if idx >= 0:
                        str = str[:idx - 1]
                    idx = str.find('-')
                    if idx >= 0:
                        str = str[:idx - 1]
                    destination_file.write(str)
                    if not str.endswith('\n'):
                        destination_file.write('\n')

        self._running = False
        self._docking_pending = False
        self._bond_pending = True

    # Add Bonds

    def _start_bonds(self):
        args = ['obabel', '-ipdb', self._ligand_output, '-osdf', '-O' + self._bond_output]
        
        nanome.util.Logs.debug("Start Bonds")
        self._start_timer = timer()
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

        docked_ligands = nanome.structure.Complex.io.from_sdf(self._bond_output)
        nanome.util.Logs.debug("Read SDF", self._bond_output)

        docked_ligands.molecular.name = "Docking"
        docked_ligands.rendering.visible = True
        if self._align == True:
            docked_ligands.transform.position = self._receptor.transform.position
            docked_ligands.transform.rotation = self._receptor.transform.rotation

        nanome.util.Logs.debug("Update workspace")
        self._plugin.make_plugin_usable()
        self._plugin.add_result_to_workspace([docked_ligands])