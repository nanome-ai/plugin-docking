import nanome
import os
import subprocess
import tempfile

from nanome.api.structure import Complex
from nanome._internal._structure._io._pdb.save import Options as _PDBOptions
from nanome.util import ComplexUtils

pdb_options = _PDBOptions()
pdb_options.write_bonds = True


class DockingCalculations():

    def __init__(self, plugin):
        self._plugin = plugin
        self.requires_site = False

    def start_docking(self, receptor, ligands, site, **params):
        docked_ligands = None
        modes = params.get('modes')
        exhaustiveness = params.get('exhaustiveness')
        align = params.get('align')

        with tempfile.TemporaryDirectory() as self.temp_dir:
            combined_ligands = ComplexUtils.combine_ligands(receptor, ligands)

            # Save all input files
            receptor_file_pdb = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb", dir=self.temp_dir)
            ligands_file_pdb = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb", dir=self.temp_dir)
            receptor.io.to_pdb(receptor_file_pdb.name, pdb_options)
            combined_ligands.io.to_pdb(ligands_file_pdb.name, pdb_options)

            # Start Ligand/ Receptor prep
            receptor_file_pdbqt = self._prepare_receptor(receptor_file_pdb)
            ligands_file_pdbqt = self._prepare_ligands(ligands_file_pdb)

            # Prepare Grid and Docking parameters.
            autogrid_input_gpf = self._prepare_grid_params(receptor_file_pdbqt, ligands_file_pdbqt, site)
            # autodock_input_dpf = self._prepare_docking_params(receptor_file_pdbqt, ligands_file_pdbqt)

            # Creates .map files and saves in the temp folder.
            self._start_autogrid4(autogrid_input_gpf)

            dock_results_pdbqt = self._start_vina(receptor_file_pdbqt, ligands_file_pdbqt, num_modes=modes, exhaustiveness=exhaustiveness)
            dock_results_sdf = self.convert_pdbqt_to_sdf(dock_results_pdbqt)
            docked_ligands = Complex.io.from_sdf(path=dock_results_sdf.name)
        
        ComplexUtils.convert_to_frames([docked_ligands])

        # make ligands invisible
        self.make_complexes_invisible(ligands)

        if len(combined_ligands.names) == 1:
            docked_ligands.name = combined_ligands.names[0] + " (Docked)"
        else:
            docked_ligands.name == "Docking Results"

        docked_ligands.visible = True
        if align:
            docked_ligands.position = receptor.position
            docked_ligands.rotation = receptor.rotation

        nanome.util.Logs.debug("Update workspace")
        self._plugin.add_result_to_workspace([docked_ligands], align)

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
        prepare_gpf4_script = os.path.join(os.path.dirname(__file__), 'prepare_gpf4.py')
        autogrid_input_gpf = tempfile.NamedTemporaryFile(delete=False, suffix=".gpf", dir=self.temp_dir)

        # gridcenter = ','.join(str(coord) for coord in site.unpack())
        grid_args = [
            'conda', 'run', '-n', 'adfr-suite',
            'python', prepare_gpf4_script,
            '-l', ligands_file_pdbqt.name,
            '-r', receptor_file_pdbqt.name,
            '-o', autogrid_input_gpf.name,
            # '-p', f"gridcenter='{gridcenter}'",
            '-y'
        ]
        subprocess.run(grid_args, cwd=self.temp_dir)
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

    def _start_vina(self, receptor_file_pdbqt, ligands_file_pdbqt, num_modes=5, exhaustiveness=8):
        # Start VINA Docking, using the autodock4 scoring.
        vina_binary = os.path.join(os.path.dirname(__file__), 'vina_1.2.2_linux_x86_64')
        dock_results = tempfile.NamedTemporaryFile(delete=False, dir=self.temp_dir, suffix='.pdbqt')
        # map files created by autogrid call, and are found using the receptor file name.
        maps_identifier = receptor_file_pdbqt.name.split('.pdbqt')[0]
        args = [
            vina_binary,
            '--scoring', 'ad4',
            '--maps', maps_identifier,
            '--ligand', ligands_file_pdbqt.name,
            '--out', dock_results.name,
            '--exhaustiveness', str(exhaustiveness),
            '--num_modes', str(num_modes)
        ]
        nanome.util.Logs.debug("Start Autodock")
        subprocess.run(args, cwd=self.temp_dir)
        return dock_results

    def make_complexes_invisible(self, complexes):
        for comp in complexes:
            comp.visible = False
        self._plugin.update_structures_shallow(complexes)

    def convert_pdbqt_to_sdf(self, pdbqt_file):
        output_file = tempfile.NamedTemporaryFile(delete=False, dir=self.temp_dir, suffix=".sdf")
        cmd = ['obabel', '-ipdbqt', pdbqt_file.name, f'-O{output_file.name}']
        subprocess.run(cmd, cwd=self.temp_dir)
        return output_file
