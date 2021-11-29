import nanome
import os
import subprocess
import sys
import tempfile

from nanome.api.structure import Complex
from nanome_docking.utils import get_complex_center


class DockingCalculations():

    def __init__(self, plugin):
        self._plugin = plugin
        self.requires_site = False

    async def start_docking(self, receptor_pdb, ligand_pdbs, site_pdb, temp_dir, **params):
        self.temp_dir = temp_dir
        modes = params.get('modes')
        exhaustiveness = params.get('exhaustiveness')
        deterministic = params.get('deterministic')

        # Get site center vector from site_pdb
        site_comp = Complex.io.from_pdb(path=site_pdb.name)
        site_center = get_complex_center(site_comp)

        # Start Ligand/ Receptor prep
        receptor_file_pdbqt = self._prepare_receptor(receptor_pdb)

        ligand_files_pdbqt = []
        for lig_pdb in ligand_pdbs:
            lig_file_pdbqt = self._prepare_ligands(lig_pdb)
            ligand_files_pdbqt.append(lig_file_pdbqt)

        output_files = []
        # Run vina, and convert output from pdbqt into a Complex object.
        for lig_file in ligand_files_pdbqt:
            # Prepare Grid and Docking parameters.
            autogrid_input_gpf = self._prepare_grid_params(receptor_file_pdbqt, lig_file, site_center)
            # autodock_input_dpf = self._prepare_docking_params(receptor_file_pdbqt, ligands_file_pdbqt)

            # Run autogrid which creates .map files and saves in the temp folder.
            self._start_autogrid4(autogrid_input_gpf)

            result_pdbqt = self._start_vina(
                receptor_file_pdbqt, lig_file, num_modes=modes, exhaustiveness=exhaustiveness,
                deterministic=deterministic)
            with open(result_pdbqt.name) as f:
                result_sdf = self.convert_pdbqt_to_sdf(f)
                output_files.append(result_sdf)
        return output_files

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
        return ligands_file_pdbqt

    def _prepare_grid_params(self, receptor_file_pdbqt, ligands_file_pdbqt, site_center):
        prepare_gpf4_script = os.path.join(os.path.dirname(__file__), 'py2', 'prepare_gpf4.py')
        autogrid_output_gpf = tempfile.NamedTemporaryFile(delete=False, suffix=".gpf", dir=self.temp_dir)

        # Write reference gpf file to set gridcenter to site
        gridcenter_line = f"gridcenter {' '.join([str(round(coord, 3)) for coord in site_center.unpack()])}"
        reference_file = tempfile.NamedTemporaryFile(suffix=".gpf", dir=self.temp_dir)
        with open(reference_file.name, 'w') as f:
            f.write(gridcenter_line)

        grid_args = [
            'conda', 'run', '-n', 'adfr-suite',
            'python', prepare_gpf4_script,
            '-l', ligands_file_pdbqt.name,
            '-r', receptor_file_pdbqt.name,
            '-o', autogrid_output_gpf.name,
            '-i', reference_file.name
        ]
        subprocess.run(grid_args, cwd=self.temp_dir)
        return autogrid_output_gpf

    def _prepare_docking_params(self, receptor_file_pdbqt, ligands_file_pdbqt):
        # Prepare Docking parameters
        prepare_dpf42_script = os.path.join(os.path.dirname(__file__), 'py2', 'prepare_dpf42.py')
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

    def _start_vina(self, receptor_file_pdbqt, ligand_file_pdbqt, num_modes=5, exhaustiveness=8, deterministic=False):
        # Start VINA Docking, using the autodock4 scoring.
        vina_binary = os.path.join(os.path.dirname(__file__), 'vina_1.2.2_linux_x86_64')
        # map files created by autogrid call, and are found using the receptor file name.
        maps_identifier = receptor_file_pdbqt.name.split('.pdbqt')[0]
        dock_results = tempfile.NamedTemporaryFile(delete=False, dir=self.temp_dir, suffix='.pdbqt')
        args = [
            vina_binary,
            '--scoring', 'ad4',
            '--maps', maps_identifier,
            '--ligand', ligand_file_pdbqt.name,
            '--out', dock_results.name,
            '--exhaustiveness', str(exhaustiveness),
            '--num_modes', str(num_modes)
        ]
        if deterministic:
            seed = '12345'
            args.extend(['--seed', seed])








        nanome.util.Logs.debug("Start Autodock")
        process = subprocess.Popen(args, cwd=self.temp_dir, stdout=subprocess.PIPE)
        self.handle_loading_bar(process, 1)
        return dock_results

    def handle_loading_bar(self, process, ligand_count):
        """Render loading bar from stdout on the menu.

        stdout has a loading bar of asterisks. Every asterisk represents about 2% completed
        """
        star_count = 0
        stars_per_complex = 51
        total_stars = stars_per_complex * ligand_count

        for c in iter(lambda: process.stdout.read(1), b''):
            if c.decode() == '*':
                star_count += 1
                self._plugin.update_loading_bar(star_count, total_stars)
            sys.stdout.buffer.write(c)

    def convert_pdbqt_to_sdf(self, pdbqt_file):
        output_file = tempfile.NamedTemporaryFile(delete=False, dir=self.temp_dir, suffix=".sdf")
        cmd = ['obabel', '-ipdbqt', pdbqt_file.name, f'-O{output_file.name}']
        subprocess.run(cmd, cwd=self.temp_dir)
        return output_file
