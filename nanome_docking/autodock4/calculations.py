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
        self.requires_site = False
        self.temp_dir = tempfile.TemporaryDirectory()


    def start_docking(self, receptor, ligands, site, exhaustiveness, modes, align, replace, scoring, visual_scores, autobox):
        temp_dir = self.temp_dir.name

        self._receptor = receptor
        self._ligands = ligands
        self._combined_ligands = ComplexUtils.combine_ligands(receptor, ligands)
        self._site = site
        self._align = align
        self._replace = replace
        self._visual_scores = visual_scores

        # Save all input files
        receptor_file_pdb = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb", dir=temp_dir)
        ligands_file_pdb = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb", dir=temp_dir)
        receptor.io.to_pdb(receptor_file_pdb.name, pdb_options)
        self._combined_ligands.io.to_pdb(ligands_file_pdb.name, pdb_options)
        assert open(receptor_file_pdb.name).read()
        assert open(receptor_file_pdb.name).read()

        # Start Ligand/ Receptor prep
        receptor_file_pdbqt = tempfile.NamedTemporaryFile(delete=False, suffix=".pdbqt", dir=temp_dir)
        ligands_file_pdbqt = tempfile.NamedTemporaryFile(delete=False, suffix=".pdbqt", dir=temp_dir)
        rec_args = [
            'conda', 'run', '-n', 'adfr-suite',
            'prepare_receptor',
            '-r', receptor_file_pdb.name,
            '-o', receptor_file_pdbqt.name,
        ]
        subprocess.run(rec_args, cwd=temp_dir)
        assert open(receptor_file_pdbqt.name).read()

        lig_args = [
            'conda', 'run', '-n', 'adfr-suite',
            'prepare_ligand',
            '-l', ligands_file_pdb.name,
            '-o', ligands_file_pdbqt.name,
            '-A', 'hydrogens',
            '-v'
        ]
        subprocess.run(lig_args, cwd=temp_dir)
        assert open(ligands_file_pdbqt.name).read()

        # Prepare Grid parameters.
        prepare_gpf4_script = os.path.join(os.path.dirname(__file__), 'prepare_gpf4.py')
        autogrid_input_gpf = tempfile.NamedTemporaryFile(delete=False, suffix=".gpf", dir=temp_dir)
        grid_args = [
            'conda', 'run', '-n', 'adfr-suite',
            'python', prepare_gpf4_script,
            '-l', ligands_file_pdbqt.name,
            '-r', receptor_file_pdbqt.name,
            '-o', autogrid_input_gpf.name
        ]
        subprocess.run(grid_args, cwd=temp_dir)
        assert open(autogrid_input_gpf.name).read()

        # Prepare Docking parameters
        prepare_dpf42_script = os.path.join(os.path.dirname(__file__), 'prepare_dpf42.py')
        autodock_input_dpf = tempfile.NamedTemporaryFile(delete=False, suffix=".dpf", dir=temp_dir)
        dock_args = [
            'conda', 'run', '-n', 'adfr-suite',
            'python', prepare_dpf42_script,
            '-l', ligands_file_pdbqt.name,
            '-r', receptor_file_pdbqt.name,
            '-o', autodock_input_dpf.name
        ]
        nanome.util.Logs.debug("Prepare grid and docking parameter files")
        subprocess.run(dock_args, cwd=temp_dir)
        assert open(autodock_input_dpf.name).read()

        # Why?
        # with open(autodock_input_dpf.name, 'a') as f:
        #     f.write("write_all\n")
        #     f.close()

        # Start Grid
        autogrid_log = tempfile.NamedTemporaryFile(delete=False, suffix=".glg", dir=temp_dir)
        args = [
            'conda', 'run', '-n', 'adfr-suite',
            'autogrid4', '-p', autogrid_input_gpf.name, '-l', autogrid_log.name
        ]
        nanome.util.Logs.debug("Start Autogrid")
        subprocess.run(args, cwd=temp_dir)

        # Start VINA Docking
        vina_binary = os.path.join(os.path.dirname(__file__), 'vina_1.2.2_linux_x86_64')
        dock_results = tempfile.NamedTemporaryFile(dir=temp_dir, suffix='.pdbqt')
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
        process = subprocess.run(args, cwd=temp_dir)

        # Start Bonds
        nanobabel_output = tempfile.NamedTemporaryFile(dir=temp_dir, suffix=".sdf")
        nanome.util.Logs.debug("Start Bonds")
        cmd = f'nanobabel convert -i {dock_results.name} -o {nanobabel_output.name}'
        args = shlex.split(cmd)
        subprocess.run(args, cwd=temp_dir)

        # make ligands invisible
        self.make_ligands_invisible()

        docked_ligands = nanome.structure.Complex.io.from_sdf(path=nanobabel_output.name)
        nanome.util.Logs.debug("Read SDF", nanobabel_output.name)
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

    def make_ligands_invisible(self):
        for ligand in self._ligands:
            ligand.visible = False
        self._plugin.update_structures_shallow(self._ligands)
