import nanome
import gzip
import itertools
import operator
import os
import subprocess
import tempfile
from timeit import default_timer as timer

# TMP
from nanome._internal._structure._io._pdb.save import Options as PDBOptions
from nanome._internal._structure._io._sdf.save import Options as SDFOptions

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
        self._pdb_options = PDBOptions()
        self._pdb_options.write_bonds = True
        self._sdf_options = SDFOptions()
        self._sdf_options.write_bonds = True

    def initialize(self):
        self._protein_input = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb")
        self._ligands_input = tempfile.NamedTemporaryFile(delete=False, suffix=".sdf")
        self._site_input = tempfile.NamedTemporaryFile(delete=False, suffix=".sdf")
        self._docking_output = tempfile.NamedTemporaryFile(delete=False, suffix=".sdf.gz")
        self._ligand_output = tempfile.NamedTemporaryFile(delete=False, suffix=".sdf")
        self._log_file = tempfile.NamedTemporaryFile(delete=False)

    def start_docking(self, receptor, ligands, site, exhaustiveness, modes, align, replace, scoring, autobox):
        self.initialize()
        
        # Save all input files
        receptor.io.to_pdb(self._protein_input.name, self._pdb_options)
        nanome.util.Logs.debug("Saved PDB", self._protein_input.name)
        ligands.io.to_sdf(self._ligands_input.name, self._sdf_options)
        nanome.util.Logs.debug("Saved SDF", self._ligands_input.name)
        site.io.to_sdf(self._site_input.name, self._sdf_options)
        nanome.util.Logs.debug("Saved SDF", self._site_input.name)

        self._exhaustiveness = exhaustiveness
        self._modes = modes
        self._receptor = receptor
        self._align = align
        self._replace = replace
        self._scoring = scoring
        self._autobox = autobox

        # Start docking process
        self._running = False
        self._request_pending = True

    def update(self):
        if self._request_pending == False:
            return

        if self._running == False:
            self._start_docking()
        elif self._check_docking():
            self._docking_finished()

    def _start_docking(self):
        if self._scoring:
            smina_args = ['./smina', '--autobox_ligand', self._site_input.name, '--score_only', '-r', self._protein_input.name, '--ligand', self._ligands_input.name, '--out', self._ligand_output.name]
        else:
            smina_args = ['./smina', '--autobox_ligand', self._site_input.name, '-r', self._protein_input.name, '--ligand', self._ligands_input.name, '--out', \
                self._docking_output.name, '--log', self._log_file.name, '--exhaustiveness', str(self._exhaustiveness), '--num_modes', str(self._modes), '--autobox_add', '+' + str(self._autobox)]
        
        nanome.util.Logs.debug("Run SMINA")
        self._start_timer = timer()
        try:
            self._smina_process = subprocess.Popen(smina_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        except:
            nanome.util.Logs.error("Couldn't execute smina, please check if executable is in the plugin folder and has permissions")
            self._request_pending = False
            self._running = False
            self._plugin.make_plugin_usable()
            self._plugin.send_notification(nanome.util.NotificationTypes.error, "Docking error, check plugin")
            return

        self._running = True
        self._plugin.send_notification(nanome.util.NotificationTypes.message, "Docking started")

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
                for result in results:
                    for line in result.split('\n'):
                        nanome.util.Logs.debug(line)
            else:
                for line in errors.splitlines():
                    nanome.util.Logs.warning(line.decode("utf-8"))
        except:
            pass

        if not self._scoring:
            self._reassemble_ligs()
        docked_ligands = nanome.structure.Complex.io.from_sdf(self._ligand_output.name)
        nanome.util.Logs.debug("Read SDF", self._ligand_output.name)

        if not self._scoring:
            docked_ligands.molecular.name = "Docking"
            docked_ligands.rendering.visible = True
            if self._align == True:
                docked_ligands.transform.position = self._receptor.transform.position
                docked_ligands.transform.rotation = self._receptor.transform.rotation
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
        os.remove(self._protein_input.name)
        os.remove(self._ligands_input.name)
        os.remove(self._site_input.name)
        os.remove(self._docking_output.name)
        os.remove(self._ligand_output.name)
        os.remove(self._log_file.name)

        self._plugin.send_notification(nanome.util.NotificationTypes.success, "Docking finished")
