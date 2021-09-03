import os
import nanome
import shlex
import subprocess
import tempfile
import re

from nanome.api.structure import Complex

PDBOPTIONS = Complex.io.PDBSaveOptions
SDFOPTIONS = Complex.io.SDFSaveOptions
from nanome.util import ComplexUtils

pdb_options = PDBOPTIONS()
pdb_options.write_bonds = True
sdf_options = SDFOPTIONS()
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

    def start_docking(self, receptor, ligands, site, exhaustiveness, modes, align, replace, scoring, visual_scores, autobox):
        self._site = site
        self._align = align
        self._replace = replace
        self._visual_scores = visual_scores
        self._receptor = receptor
        self._ligands = ligands

        docked_ligands = None
        with tempfile.TemporaryDirectory() as temp_dir:
            self.temp_dir = temp_dir
            combined_ligands = ComplexUtils.combine_ligands(receptor, ligands)

            # Save all input files
            receptor_file_pdb = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb", dir=temp_dir)
            ligands_file_pdb = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb", dir=temp_dir)
            receptor.io.to_pdb(receptor_file_pdb.name, pdb_options)
            combined_ligands.io.to_pdb(ligands_file_pdb.name, pdb_options)
            assert open(receptor_file_pdb.name).read()
            assert open(receptor_file_pdb.name).read()

            # Start Ligand/ Receptor prep
            receptor_file_pdbqt = self._prepare_receptor(receptor_file_pdb)
            ligands_file_pdbqt = self._prepare_ligands(ligands_file_pdb)
            
            # Prepare Grid parameters.
            autogrid_input_gpf = self._prepare_grid_params(receptor_file_pdbqt, ligands_file_pdbqt)            
            autodock_input_dpf = self._prepare_docking_params(receptor_file_pdbqt, ligands_file_pdbqt)

            grid_filepaths = self._start_autogrid4(autogrid_input_gpf)

            dock_results = self._start_vina(receptor_file_pdbqt, ligands_file_pdbqt)

            docked_ligands_sdf = self._calculate_bonds(dock_results)
            results_complex = nanome.structure.Complex.io.from_sdf(path=docked_ligands_sdf.name)
        
        # make ligands invisible
        self.make_ligands_invisible()

        if len(combined_ligands.names) == 1:
            results_complex.name = combined_ligands.names[0] + " (Docked)"
        else:
            results_complex.name == "Docking Results"

        results_complex.visible = True
        if self._align is True:
            results_complex.transform.position = self._receptor.transform.position
            results_complex.transform.rotation = self._receptor.transform.rotation

        nanome.util.Logs.debug("Update workspace")
        # TODO: verify this shouldn't be here anymore (test)
        # self._plugin.make_plugin_usable()
        self._plugin.add_result_to_workspace([results_complex])

    def _prepare_receptor(self, pdb_file):
        """Convert pdb file into pdbqt."""
        receptor_file_pdbqt = tempfile.NamedTemporaryFile(delete=False, suffix=".pdbqt", dir=self.temp_dir)
        rec_args = [
            'conda', 'run', '-n', 'adfr-suite',
            'prepare_receptor',
            '-r', pdb_file.name,
            '-o', receptor_file_pdbqt.name,
        ]
        subprocess.run(rec_args, cwd=self.temp_dir)
        assert open(receptor_file_pdbqt.name).read()
        return receptor_file_pdbqt
    
    def _prepare_ligands(self, ligands_file_pdb):
        """Convert pdb file into pdbqt."""
        ligands_file_pdbqt = tempfile.NamedTemporaryFile(delete=False, suffix=".pdbqt", dir=self.temp_dir)
        lig_args = [
            'conda', 'run', '-n', 'adfr-suite',
            'prepare_ligand',
            '-l', ligands_file_pdb.name,
            '-o', ligands_file_pdbqt.name,
            '-A', 'hydrogens',
            '-v'
        ]
        subprocess.run(lig_args, cwd=self.temp_dir)
        assert open(ligands_file_pdbqt.name).read()
        return ligands_file_pdbqt

    def _prepare_grid_params(self, receptor_file_pdbqt, ligands_file_pdbqt):
        prepare_gpf4_script = os.path.join(os.path.dirname(__file__), 'prepare_gpf4.py')
        autogrid_input_gpf = tempfile.NamedTemporaryFile(delete=False, suffix=".gpf", dir=self.temp_dir)
        grid_args = [
            'conda', 'run', '-n', 'adfr-suite',
            'python', prepare_gpf4_script,
            '-l', ligands_file_pdbqt.name,
            '-r', receptor_file_pdbqt.name,
            '-o', autogrid_input_gpf.name
        ]
        subprocess.run(grid_args, cwd=self.temp_dir)
        assert open(autogrid_input_gpf.name).read()
        return autogrid_input_gpf

    def _prepare_docking_params(self, receptor_file_pdbqt, ligands_file_pdbqt):
        # Prepare Docking parameters
        prepare_dpf42_script = os.path.join(os.path.dirname(__file__), 'prepare_dpf42.py')
        autodock_input_dpf = tempfile.NamedTemporaryFile(delete=False, suffix=".dpf", dir=self.temp_dir)
        dock_args = [
            'conda', 'run', '-n', 'adfr-suite',
            'python', prepare_dpf42_script,
            '-l', ligands_file_pdbqt.name,
            '-r', receptor_file_pdbqt.name,
            '-o', autodock_input_dpf.name
        ]
        nanome.util.Logs.debug("Prepare grid and docking parameter files")
        subprocess.run(dock_args, cwd=self.temp_dir)
        assert open(autodock_input_dpf.name).read()
        return autodock_input_dpf

    def _start_autogrid4(self, autogrid_input_gpf):
        # Start Grid
        autogrid_log = tempfile.NamedTemporaryFile(delete=False, suffix=".glg", dir=self.temp_dir)
        args = [
            'conda', 'run', '-n', 'adfr-suite',
            'autogrid4', '-p', autogrid_input_gpf.name, '-l', autogrid_log.name
        ]
        nanome.util.Logs.debug("Start Autogrid")
        subprocess.run(args, cwd=self.temp_dir)
        generated_filepaths = [
            f'{self.temp_dir}/{filename}' for filename in os.listdir(self.temp_dir)
            if filename.endswith('.map') or filename.endswith('.fld')
        ]
        return generated_filepaths

    def _start_vina(self, receptor_file_pdbqt, ligands_file_pdbqt):
        # Start VINA Docking
        vina_binary = os.path.join(os.path.dirname(__file__), 'vina_1.2.2_linux_x86_64')
        dock_results = tempfile.NamedTemporaryFile(dir=self.temp_dir, suffix='.pdbqt')
        args = [
            vina_binary,
            '--receptor', receptor_file_pdbqt.name,
            '--ligand', ligands_file_pdbqt.name,
            '--center_x', '3.37',
            '--center_y', '-17.2',
            '--center_z', '8.49',
            '--size_x', '4',
            '--size_y', '4',
            '--size_z', '4',
            '--out', dock_results.name
        ]
        nanome.util.Logs.debug("Start Autodock")
        process = subprocess.run(args, cwd=self.temp_dir)
        assert open(dock_results.name).read()
        return dock_results

    def _calculate_bonds(self, dock_results):
        # Start Bonds
        nanobabel_output = tempfile.NamedTemporaryFile(dir=self.temp_dir, suffix=".sdf")
        nanome.util.Logs.debug("Start Bonds")
        cmd = f'nanobabel convert -i {dock_results.name} -o {nanobabel_output.name}'
        args = shlex.split(cmd)
        subprocess.run(args, cwd=self.temp_dir)
        assert open(nanobabel_output.name).read()
        return nanobabel_output

    def make_ligands_invisible(self):
        for ligand in self._ligands:
            ligand.visible = False
        self._plugin.update_structures_shallow(self._ligands)
