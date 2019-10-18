import traceback

import nanome
import gzip
import itertools
import operator
import os
import subprocess
import tempfile
import stat
from timeit import default_timer as timer

from nanome.util.enums import NotificationTypes

SDFOPTIONS = nanome.api.structure.Complex.io.SDFSaveOptions()
SDFOPTIONS.write_bonds = True
PDBOPTIONS = nanome.api.structure.Complex.io.PDBSaveOptions()
PDBOPTIONS.write_bonds = True

SMINA_PATH = os.path.join(os.path.dirname(__file__), 'smina')
try:
    os.chmod(SMINA_PATH, stat.S_IXGRP | stat.S_IEXEC)
except:
    pass

def not_dollars(line):
    return b'$$$$' != line.strip(b'\n')

def parse_molecules(file):
    nanome.util.Logs.debug("Parsing output", file)
    with gzip.open(file) as lines:
        while True:
            block = list(itertools.takewhile(not_dollars, lines ))
            if not block:
                break
            yield block + [b'$$$$\n']

class DockingCalculations():
    def __init__(self, plugin):
        self._plugin = plugin
        self._request_pending = False
        self._running = False
        self._started_conversion = False
        self._started_docking = False

    def initialize(self):
        self._receptor_input = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb")
        self._ligands_input = tempfile.NamedTemporaryFile(delete=False, suffix=".sdf")
        self._ligands_pdb_input = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb")
        self._site_input = tempfile.NamedTemporaryFile(delete=False, suffix=".sdf")
        self._docking_output = tempfile.NamedTemporaryFile(delete=False, suffix=".sdf.gz")
        self._ligand_output = tempfile.NamedTemporaryFile(delete=False, suffix=".sdf")
        self._log_file = tempfile.NamedTemporaryFile(delete=False)

    def start_docking(self, receptor, ligands, site, exhaustiveness, modes, align, replace, scoring, autobox):
        self._exhaustiveness = exhaustiveness
        self._modes = modes
        self._receptor = receptor
        self._ligands = ligands
        self._site = site
        self._align = align
        self._replace = replace
        self._scoring = scoring
        self._autobox = autobox

        # Start docking process
        self._running = False
        self._started_conversion = False
        self._request_pending = True

    def update(self):
        if self._request_pending == False:
            return

        if not self._running and not self._started_conversion:
            self._start_converting()
            self._started_conversion = True
        elif self._check_conversion() and not self._started_docking:
            self._start_docking()
            self._started_docking = True
        elif self._check_docking():
            self._docking_finished()

    def _start_converting(self):
        self.initialize()
         # Save all input files
        self._receptor.io.to_pdb(self._receptor_input.name, PDBOPTIONS)
        nanome.util.Logs.debug("Saved PDB", self._receptor_input.name)
        self._ligands.io.to_sdf(self._ligands_input.name, SDFOPTIONS)
        nanome.util.Logs.debug("Saved SDF", self._ligands_input.name)
        self._site.io.to_sdf(self._site_input.name, SDFOPTIONS)
        nanome.util.Logs.debug("Saved SDF", self._site_input.name)

        args = ['obabel', '-isdf', self._ligands_input.name, '-opdb', '-O' + self._ligands_pdb_input.name]
        self._obabel_process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
    def _check_conversion(self):
        poll = self._obabel_process.poll()
        if poll == None:
            self._obabel_process.communicate()
        return poll != None

    def _start_docking(self):
        exe_path = os.path.join(os.path.dirname(__file__), 'smina')
        if self._scoring:
            smina_args = [exe_path, '--autobox_ligand', self._site_input.name, '--score_only', '-r', self._receptor_input.name, '--ligand', self._ligands_input.name, '--out', self._ligand_output.name]
        else:
            smina_args = [exe_path, '--autobox_ligand', self._site_input.name, '-r', self._receptor_input.name, '--ligand', self._ligands_pdb_input.name, '--out', \
                self._docking_output.name, '--log', self._log_file.name, '--exhaustiveness', str(self._exhaustiveness), '--num_modes', str(self._modes), '--autobox_add', '+' + str(self._autobox)]

        nanome.util.Logs.debug("Run SMINA")
        self._start_timer = timer()
        try:
            self._smina_process = subprocess.Popen(smina_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        except:
            nanome.util.Logs.error("Couldn't execute smina, please check if executable is in the plugin folder and has permissions. Path:", SMINA_PATH)
            self._request_pending = False
            self._running = False
            self._plugin.make_plugin_usable()
            self._plugin.send_notification(NotificationTypes.error, "Docking error, check plugin")
            return

        self._running = True
        self._plugin.send_notification(NotificationTypes.message, "Docking started")

    def _check_docking(self):
        return self._smina_process.poll() != None
                        
    def _reassemble_ligs(self):
        generators = parse_molecules(self._docking_output.name)

        fout = open(self._ligand_output.name, 'w')
            
        try:
            for gen in generators:
                for line in gen:
                    str = line.decode("utf-8").strip('\n')
                    fout.write(''.join(str))
                    fout.write('\n')
        except StopIteration:
            pass
                
        fout.close()
    
    def _docking_finished(self):
        end = timer()
        nanome.util.Logs.debug("Docking Finished in", end - self._start_timer, "seconds")
        self._request_pending = False
        try:
            (results, errors) = self._smina_process.communicate()
            if len(errors) == 0:
                for line in results.decode().splitlines():
                    nanome.util.Logs.debug(line)
            else:
                for line in errors.decode().splitlines():
                    nanome.util.Logs.warning(line.decode("utf-8"))
        except Exception as e:
            print(traceback.format_exc())

        if not self._scoring:
            self._reassemble_ligs()
        docked_ligands = nanome.structure.Complex.io.from_sdf(path=self._ligand_output.name)
        nanome.util.Logs.debug("Read SDF", self._ligand_output.name)

        if not self._scoring:
            docked_ligands.molecular.name = "Docking"
            docked_ligands.rendering.visible = True
            if self._align == True:
                mat = self._site.transform.get_complex_to_workspace_matrix()
                docked_ligands.transform.position = self._site.transform.position
                docked_ligands.transform.rotation = self._site.transform.rotation
                self._site.transform.position = self._site.transform.position
                docked_ligands.rendering.boxed = True

        self._plugin.make_plugin_usable()
        if self._scoring:
            nanome.util.Logs.debug("Display scoring result")
            self._plugin.display_scoring_result(docked_ligands)
        else:
            nanome.util.Logs.debug("Update workspace")
            self._plugin.add_result_to_workspace([docked_ligands])

        self._docking_output.close()
        self._ligand_output.close()
        self._log_file.close()
        os.remove(self._receptor_input.name)
        os.remove(self._ligands_input.name)
        os.remove(self._site_input.name)
        os.remove(self._docking_output.name)
        os.remove(self._ligand_output.name)
        os.remove(self._log_file.name)

        self._plugin.send_notification(NotificationTypes.success, "Docking finished")
