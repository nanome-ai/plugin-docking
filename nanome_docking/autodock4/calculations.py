import os
import nanome
import shlex
import subprocess
import tempfile
import re

# TMP
from nanome._internal._structure._io._pdb.save import Options as PDBOptions
from nanome._internal._structure._io._sdf.save import Options as SDFOptions
from nanome.util import ComplexUtils

pdb_options = PDBOptions()
pdb_options.write_bonds = True
sdf_options = SDFOptions()
sdf_options.write_bonds = True


class DockingCalculations():

    def __init__(self, plugin):
        self._plugin = plugin
        self._preparation_pending = False
        self._grid_pending = False
        self._docking_pending = False
        self._bond_pending = False
        self._running = False
        self.requires_site = False

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
        self.requires_site = False

    def start_docking(self, receptor, ligands, site, exhaustiveness, modes, align, replace, scoring, visual_scores, autobox):
        self._receptor = receptor
        self._ligands = ligands
        self._combined_ligands = ComplexUtils.combine_ligands(receptor, ligands)
        self._site = site
        self._align = align
        self._replace = replace
        self._visual_scores = visual_scores

        # Save all input files
        receptor.io.to_pdb(self._protein_input.name, pdb_options)
        nanome.util.Logs.debug("Saved PDB", self._protein_input.name)
        self._combined_ligands.io.to_pdb(self._ligands_input.name, pdb_options)
        nanome.util.Logs.debug("Saved PDB", self._ligands_input.name)

        # Start docking process
        self._prepare_receptor(self._protein_input.name, self._protein_input_converted.name)
        self._prepare_ligands(self._ligands_input.name, self._ligands_input_converted.name)

        # Prepare Grid and Docking parameters.
        print("protein input converted name:", self._protein_input_converted.name)
        grid_args = [
            'conda', 'run', '-n', 'adfr-suite',
            'python', os.path.join(os.path.dirname(__file__), 'prepare_gpf4.py'),
            '-l', self._ligands_input_converted.name,
            '-r', self._protein_input_converted.name,
            '-o', self._autogrid_input.name
        ]
        dock_args = [
            'conda', 'run', '-n', 'adfr-suite',
            'python', os.path.join(os.path.dirname(__file__), 'prepare_dpf42.py'),
            '-l', self._ligands_input_converted.name,
            '-r', self._protein_input_converted.name,
            '-o', self._autodock_input.name
        ]

        nanome.util.Logs.debug("Prepare grid and docking parameter files")
        subprocess.run(grid_args, cwd=self.temp_dir.name)
        subprocess.run(dock_args, cwd=self.temp_dir.name)
        assert open(self._autodock_input.name).read()
        assert open(self._autogrid_input.name).read()

        # Why?
        with open(self._autodock_input.name, 'a') as f:
            f.write("write_all\n")
            f.close()

        # Start Grid
        param_filename = self._autogrid_input.name
        full_name_log = self._autogrid_log.name
        path = os.path.dirname(param_filename)
        args = [
            'conda', 'run', '-n', 'adfr-suite',
            'autogrid4', '-p', param_filename, '-l', full_name_log
        ]
        nanome.util.Logs.debug("Start Autogrid")
        subprocess.run(args, cwd=path)

        # Start Docking
        full_name_input = self._autodock_input.name
        path = os.path.dirname(full_name_input)
        full_name_log = self._autodock_log.name

        vina_binary = os.path.join(os.path.dirname(__file__), 'vina_1.2.2_linux_x86_64')
        dock_results = tempfile.NamedTemporaryFile(dir=self.temp_dir.name)
        args = [
            vina_binary,
            '--receptor', self._protein_input_converted.name,
            '--ligand', self._ligands_input_converted.name,
            '--center_x', '3.37',
            '--center_y', '-17.2',
            '--center_z', '8.49',
            '--size_x', '4',
            '--size_y', '4',
            '--size_z', '4',
            '--out', dock_results
        ]   

        nanome.util.Logs.debug("Start Autodock")
        process = subprocess.run(args, cwd=path)

        # Start Bonds
        nanome.util.Logs.debug("Start Bonds")
        cmd = f'nanobabel convert -i {dock_results} -o {self._bond_output.name}'
        args = shlex.split(cmd)
        self._nanobabel_process = subprocess.run(args, cwd=path)

        # make ligands invisible
        self.make_ligands_invisible()

        docked_ligands = nanome.structure.Complex.io.from_sdf(path=self._bond_output.name)
        nanome.util.Logs.debug("Read SDF", self._bond_output.name)
        if len(self._combined_ligands.names) == 1:
            docked_ligands.name = self._combined_ligands.names[0] + " (Docked)"
        else:
            docked_ligands.name == "Docking Results"

        docked_ligands.visible = True
        if self._align is True:
            docked_ligands.transform.position = self._receptor.transform.position
            docked_ligands.transform.rotation = self._receptor.transform.rotation

        nanome.util.Logs.debug("Update workspace")
        # TODO: verify this shouldn't be here anymore (test)
        # self._plugin.make_plugin_usable()
        self._plugin.add_result_to_workspace([docked_ligands])

    def _prepare_ligands(self, input_file, output_file):
        lig_args = [
            'conda', 'run', '-n', 'adfr-suite',
            'prepare_ligand',
            '-l', input_file,
            '-o', output_file,
            '-A', 'hydrogens',
            '-v'
        ]
        process = subprocess.run(lig_args, cwd=self.temp_dir.name)
        assert open(output_file).read()
        return process

    def _prepare_receptor(self, input_file, output_file):
        rec_args = [
            'conda', 'run', '-n', 'adfr-suite',
            'prepare_receptor',
            '-r', input_file,
            '-o', output_file,
        ]
        subprocess.run(rec_args, cwd=self.temp_dir.name)
        assert open(output_file).read()

    # Autogrid
    def _start_grid(self):
        param_filename = self._autogrid_input.name
        full_name_log = self._autogrid_log.name
        path = os.path.dirname(param_filename)

        args = [
            'conda', 'run', '-n', 'adfr-suite',
            'autogrid4', '-p', param_filename, '-l', full_name_log]

        nanome.util.Logs.debug("Start Autogrid")
        self._autogrid_process = subprocess.run(args, cwd=path)
        self._running = True

    def _check_grid(self):
        return self._autogrid_process.poll() is not None

    def make_ligands_invisible(self):
        for ligand in self._ligands:
            ligand.visible = False
        self._plugin.update_structures_shallow(self._ligands)
