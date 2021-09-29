import nanome
import os
import subprocess
import tempfile

from nanome.api.structure import Complex
from nanome_docking.utils import get_complex_center
from nanome._internal._structure._io._pdb.save import Options as _PDBOptions
from nanome.util import ComplexUtils

pdb_options = _PDBOptions()
pdb_options.write_bonds = True


class DockingCalculations():

    def __init__(self, plugin):
        self._plugin = plugin
        self.requires_site = False

    async def start_docking(self, receptor, ligands, site, **params):
        docked_ligands = []
        modes = params.get('modes')
        exhaustiveness = params.get('exhaustiveness')
        align = params.get('align')

        with tempfile.TemporaryDirectory() as self.temp_dir, tempfile.TemporaryDirectory() as self.output_dir:
            # Save all input files
            receptor_file_pdb = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb", dir=self.temp_dir)
            receptor.io.to_pdb(receptor_file_pdb.name, pdb_options)

            # Map ligand index to its pdb file, so that we can find it later
            lig_files_pdb = []
            for lig in ligands:
                lig_pdb = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb", dir=self.temp_dir)
                lig.io.to_pdb(lig_pdb.name, pdb_options)
                lig_files_pdb.append((lig.index, lig_pdb))

            # Start Ligand/ Receptor prep
            receptor_file_pdbqt = self._prepare_receptor(receptor_file_pdb)

            ligand_tracker = []
            for index, lig_pdb in lig_files_pdb:
                lig_file_pdbqt = self._prepare_ligands(lig_pdb)
                ligand_tracker.append((index, lig_file_pdbqt))

            # Prepare Grid and Docking parameters.
            autogrid_input_gpf = self._prepare_grid_params(receptor_file_pdbqt, ligand_tracker[0][1], site)
            # autodock_input_dpf = self._prepare_docking_params(receptor_file_pdbqt, ligands_file_pdbqt)

            # Creates .map files and saves in the temp folder.
            self._start_autogrid4(autogrid_input_gpf)

            # Run vina, and convert output from pdbqt into a Complex object.
            ligand_files = [f for index, f in ligand_tracker]
            dock_results_dir = self._start_vina(
                receptor_file_pdbqt, ligand_files, self.output_dir, num_modes=modes, exhaustiveness=exhaustiveness)

            for dock_result in os.listdir(self.output_dir):
                filepath = f'{dock_results_dir}/{dock_result}'
                with open(filepath) as f:
                    # Look up original ligand complex
                    result_filename = dock_result.split('_out')[0]
                    comp_index = next(
                        index for index, pdbqt_file in ligand_files if result_filename in pdbqt_file.name)
                    original_lig = next(lig for lig in ligands if lig.index == comp_index)
                    # Convert pdbqt into new Complex object.
                    dock_results_sdf = self.convert_pdbqt_to_sdf(f)
                    docked_ligand = Complex.io.from_sdf(path=dock_results_sdf.name)
                    docked_ligand.full_name = f'{original_lig.full_name} (Docked)'
                    docked_ligands.append(docked_ligand)

        ComplexUtils.convert_to_frames(docked_ligands)

        # Make original ligands hidden, and add docked ligands to workspace.
        self.make_complexes_invisible(ligands)
        for comp in docked_ligands:
            if align:
                comp.position = receptor.position
                comp.rotation = receptor.rotation
            comp.visible = True
        nanome.util.Logs.debug("Update workspace")
        self._plugin.add_result_to_workspace(docked_ligands, align)

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

    def _prepare_grid_params(self, receptor_file_pdbqt, ligands_file_pdbqt, site):
        prepare_gpf4_script = os.path.join(os.path.dirname(__file__), 'py2', 'prepare_gpf4.py')
        autogrid_output_gpf = tempfile.NamedTemporaryFile(delete=False, suffix=".gpf", dir=self.temp_dir)

        # Write reference gpf file to set gridcenter to site
        site_center = get_complex_center(site)
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

    def _start_vina(self, receptor_file_pdbqt, ligand_files_pdbqt, output_dir, num_modes=5, exhaustiveness=8):
        # Start VINA Docking, using the autodock4 scoring.
        vina_binary = os.path.join(os.path.dirname(__file__), 'vina_1.2.2_linux_x86_64')
        # map files created by autogrid call, and are found using the receptor file name.
        maps_identifier = receptor_file_pdbqt.name.split('.pdbqt')[0]

        batch_args = []
        for lig_file in ligand_files_pdbqt:
            batch_args += ['--batch', lig_file.name]
        args = [
            vina_binary,
            '--scoring', 'ad4',
            '--maps', maps_identifier,
            *batch_args,
            '--dir', output_dir,
            '--exhaustiveness', str(exhaustiveness),
            '--num_modes', str(num_modes)
        ]
        nanome.util.Logs.debug("Start Autodock")
        process = subprocess.Popen(args, cwd=self.temp_dir, stdout=subprocess.PIPE)

        # stdout has a loading bar of asterisks. Every asterisk represents about 2% completed
        # update loading bar on menu accordingly
        star_count = 0
        stars_per_complex = 51
        total_stars = stars_per_complex * len(ligand_files_pdbqt)

        for c in iter(lambda: process.stdout.read(1), b''):
            if c.decode() == '*':
                star_count += 1
                self._plugin.update_loading_bar(star_count, total_stars)
            # import sys
            # sys.stdout.buffer.write(c)
        return output_dir

    def make_complexes_invisible(self, complexes):
        for comp in complexes:
            comp.visible = False
        self._plugin.update_structures_shallow(complexes)

    def convert_pdbqt_to_sdf(self, pdbqt_file):
        output_file = tempfile.NamedTemporaryFile(delete=False, dir=self.temp_dir, suffix=".sdf")
        cmd = ['obabel', '-ipdbqt', pdbqt_file.name, f'-O{output_file.name}']
        subprocess.run(cmd, cwd=self.temp_dir)
        return output_file
