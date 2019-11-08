import traceback
import math
import re
import nanome
import shutil
import gzip
import itertools
import operator
import os
import subprocess
import tempfile
import stat
from timeit import default_timer as timer

from nanome.util.enums import NotificationTypes

from .ComplexUtils import ComplexUtils

DEBUG = True

SDFOPTIONS = nanome.api.structure.Complex.io.SDFSaveOptions()
SDFOPTIONS.write_bonds = True
PDBOPTIONS = nanome.api.structure.Complex.io.PDBSaveOptions()
PDBOPTIONS.write_bonds = True

SMINA_PATH = os.path.join(os.path.dirname(__file__), 'smina')
try:
    os.chmod(SMINA_PATH, stat.S_IXGRP | stat.S_IEXEC)
except:
    pass

class DockingCalculations():
    def __init__(self, plugin):
        self._plugin = plugin
        self._request_pending = False
        self._running = False
        self._structures_written = False
        self._started_docking = False
        self.requires_site = True


    def initialize(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self._receptor_input = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb", dir=self.temp_dir.name)
        self._ligands_input = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb", dir=self.temp_dir.name)
        self._site_input = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb", dir=self.temp_dir.name)
        self._docking_output = tempfile.NamedTemporaryFile(delete=False, prefix="output", suffix=".sdf", dir=self.temp_dir.name)
        self._ligand_output = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb", dir=self.temp_dir.name)
        self._log_file = tempfile.NamedTemporaryFile(delete=False, dir=self.temp_dir.name)

    def start_docking(self, receptor, ligands, site, exhaustiveness, modes, align, replace, scoring, visual_scores, autobox):
        self._exhaustiveness = exhaustiveness
        self._modes = modes
        self._receptor = receptor
        self._combined_ligands = ComplexUtils.combine_ligands(receptor, ligands)
        self._ligands = ligands
        self._site = site
        self._align = align
        self._replace = replace
        self._scoring = scoring
        self._visual_scores = visual_scores
        self._autobox = autobox

        # Start docking process
        self._running = False
        self._structures_written = False
        self._started_docking = False
        self._request_pending = True


    def update(self):
        if self._request_pending == False:
            return

        if not self._running and not self._structures_written:
            self.write_structures_to_file()
            self._structures_written = True
        elif not self._started_docking:
            self._start_docking()
            self._started_docking = True
        elif self._check_docking():
            self._docking_finished()

    def write_structures_to_file(self):
        self.initialize()
         # Save all input files
        self._receptor.io.to_pdb(self._receptor_input.name, PDBOPTIONS)
        nanome.util.Logs.debug("Saved PDB", self._receptor_input.name)
        self._combined_ligands.io.to_pdb(self._ligands_input.name, PDBOPTIONS)
        nanome.util.Logs.debug("Saved PDB", self._ligands_input.name)
        self._site.io.to_pdb(self._site_input.name, PDBOPTIONS)
        nanome.util.Logs.debug("Saved PDB", self._site_input.name)
        
    def _check_conversion(self):
        poll = self._obabel_process.poll()
        if poll == None:
            self._obabel_process.communicate()
        return poll != None

    def _start_docking(self):
        exe_path = os.path.join(os.path.dirname(__file__), 'smina')
        
        if self._scoring:
            smina_args = [exe_path, '-r', self._receptor_input.name, '-l', self._ligands_input.name, '--autobox_ligand', self._site_input.name, '--score_only', '--out', self._ligand_output.name]
        else:

            smina_args = [exe_path, '-r', self._receptor_input.name, '-l', self._ligands_input.name, '--autobox_ligand', self._site_input.name, '--out', \
                self._docking_output.name, '--log', self._log_file.name, '--exhaustiveness', str(self._exhaustiveness), '--num_modes', str(self._modes), '--autobox_add', str(self._autobox), '--atom_term_data']
        
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

    def _resume_docking_finished(self, docking_results):
        docking_results = docking_results[0]
        nanome.util.Logs.debug("Read SDF", self._docking_output.name)
        for i, molecule in enumerate(docking_results.molecules):
            nanome.util.Logs.debug("molecule " + str(i))
            self.set_scores(i, molecule)

        docking_results.set_current_frame(0)

        if self._visual_scores:
            self.visualize_scores(docking_results)
            
        if self._scoring:
            nanome.util.Logs.debug("Display scoring result")
            self._plugin.display_scoring_result(docking_results)
        else:
            nanome.util.Logs.debug("Update workspace")
            self._plugin.add_result_to_workspace([docking_results], self._align)

        shutil.rmtree(self.temp_dir.name)

        self._plugin.send_notification(NotificationTypes.success, "Docking finished")

    def make_ligands_invisible(self):
        for ligand in self._ligands:
            ligand.visible = False

        self._plugin.update_structures_shallow(self._ligands)
    
    def _docking_finished(self):
        end = timer()
        nanome.util.Logs.debug("Docking Finished in", end - self._start_timer, "seconds")
        self._request_pending = False

        self.make_ligands_invisible()

        docking_results = nanome.structure.Complex.io.from_sdf(path=self._docking_output.name)
        docking_results.index = -1
        if not self._scoring:
            docking_results.visible = True
            docking_results.locked = True
            if len(self._combined_ligands.names) > 1:
                docking_results.name += "Docking Results"
            elif len(self._combined_ligands.names) == 1:
                docking_results.name = self._combined_ligands.names[0] + " (Docked)"
            
        self._plugin.replace_conformer([docking_results], self._resume_docking_finished, existing=False)

    def update_min_max_scores(self, molecule, score):
        min_score = molecule.min_atom_score
        max_score = molecule.max_atom_score
        molecule.min_atom_score = score if score < min_score else min_score
        molecule.max_atom_score = score if score > max_score else max_score

    def set_scores(self, lig_i, molecule):
        molecule.min_atom_score = float('inf')
        molecule.max_atom_score = float('-inf')

        num_rgx = '(-?[\d.]+(?:e[+-]\d+)?)'
        pattern = re.compile('<{},{},{}> {} {} {} {} {}'.format(*([num_rgx] * 8)), re.U)
        for associated in molecule.associateds:
            pose_score = associated['> <minimizedAffinity>']
            for residue in molecule.residues:
                residue.label_text = pose_score + " kcal/mol"
                residue.labeled = True
            interaction_terms = associated['> <atomic_interaction_terms>']
            interaction_values = re.findall(pattern, interaction_terms)
            atom_count = len([atom for atom in molecule.atoms])
            for i, atom in enumerate(molecule.atoms):
                if i < len(interaction_values) - 1:
                    nanome.util.Logs.debug("interaction values for atom " + str(i) + ": "+ str(interaction_values[i]))
                    atom.score = float(interaction_values[i][5])
                    self.update_min_max_scores(molecule, atom.score)

    
    def visualize_scores(self, ligand_complex):
        for molecule in ligand_complex.molecules:
            for atom in molecule.atoms:
                if hasattr(atom, "score"):
                    atom.atom_mode = nanome.api.structure.Atom.AtomRenderingMode.Point
                    denominator = -molecule.min_atom_score if atom.score < 0 else molecule.max_atom_score
                    norm_score = atom.score / denominator
                    atom.atom_scale = norm_score * 1.5 + 0.1
                    atom.label_text = self.truncate(atom.score, 3)
                    # atom.labeled = True

    def truncate(self, f, n):
        '''Truncates/pads a float f to n decimal places without rounding'''
        s = '{}'.format(f)
        if 'e' in s or 'E' in s:
            return '{0:.{1}f}'.format(f, n)
        i, p, d = s.partition('.')
        return '.'.join([i, (d+'0'*n)[:n]])