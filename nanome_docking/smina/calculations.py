import sys
import time
import os
import tempfile
import subprocess
from nanome.util import Logs

SMINA_PATH = os.path.join(os.getcwd(), 'nanome_docking', 'smina', 'smina_binary')


class DockingCalculations():

    def __init__(self, plugin):
        self.plugin = plugin
        self.requires_site = True
        self.loading_bar_counter = 0

    async def start_docking(self, receptor_pdb, ligand_pdbs, site_pdb, temp_dir, exhaustiveness=None, modes=None, autobox=None, deterministic=None, **kwargs):
        # Start docking process
        start_time = time.time()
        Logs.message("Smina Calculation started.")
        self.loading_bar_counter = 0
        log_file = tempfile.NamedTemporaryFile(delete=False, dir=temp_dir)
        smina_output_sdfs = []
        for ligand_pdb in ligand_pdbs:
            # Read first line to get the number of frames
            nummdl_line = ligand_pdb.readline().decode()
            if nummdl_line.startswith("NUMMDL"):
                frame_count = int(nummdl_line.split()[1])
            else:
                Logs.warning("NUMMDL line not found in PDB file. Assuming 1 frame.")
                frame_count = 1
            output_sdf = tempfile.NamedTemporaryFile(delete=False, prefix="output", suffix=".sdf", dir=temp_dir)
            process = self.run_smina(ligand_pdb, receptor_pdb, site_pdb, output_sdf, log_file, exhaustiveness, modes, autobox, frame_count, deterministic)
            self.handle_loading_bar(process, frame_count)
            smina_output_sdfs.append(output_sdf)
        end_time = time.time()
        Logs.message("Smina Calculation finished in {} seconds.".format(round(end_time - start_time, 2)))
        process.terminate()
        return smina_output_sdfs

    def run_smina(self, ligand_pdb, receptor_pdb, site_pdb, output_sdf, log_file,
                  exhaustiveness=None, modes=None, autobox=None, ligand_count=1,
                  deterministic=False, **kwargs):
        smina_args = [
            '-r', receptor_pdb.name,
            '-l', ligand_pdb.name,
            '--autobox_ligand', site_pdb.name,
            '--out', output_sdf.name,
            '--log', log_file.name,
            '--exhaustiveness', str(exhaustiveness),
            '--num_modes', str(modes),
            '--autobox_add', str(autobox),
            '--atom_term_data'
        ]

        # To make runs deterministic, we manually set the seed. Otherwise random seed is used.
        if deterministic:
            seed = 0
            smina_args.extend(['--seed', seed])

        cmd = [SMINA_PATH, *smina_args]
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        self.handle_loading_bar(process, ligand_count)
        process.terminate()
        return process

    def handle_loading_bar(self, process, frame_count):
        """Render loading bar from stdout on the menu.

        stdout has a loading bar of asterisks. Every asterisk represents about 2% completed
        """
        stars_per_complex = 51
        total_stars = stars_per_complex * frame_count
        self.loading_bar_counter = 0
        for c in iter(lambda: process.stdout.read(1), b''):
            if c.decode() == '*':
                self.loading_bar_counter += 1
                self.plugin.update_loading_bar(self.loading_bar_counter, total_stars)
            sys.stdout.buffer.write(c)
