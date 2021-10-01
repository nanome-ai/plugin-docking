import sys
import nanome
import os
import tempfile
import subprocess
from nanome.util import ComplexUtils, Logs


PDBOPTIONS = nanome.api.structure.Complex.io.PDBSaveOptions()
PDBOPTIONS.write_bonds = True

SMINA_PATH = os.path.join(os.getcwd(), 'nanome_docking', 'smina', 'smina_binary')


class DockingCalculations():

    def __init__(self, plugin):
        self.plugin = plugin
        self.requires_site = True
        self.loading_bar_counter = 0

    async def start_docking(self, receptor, ligands, site, temp_dir,     exhaustiveness=None, modes=None, autobox=None, **kwargs):
        # Start docking process
        receptor_pdb = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb", dir=temp_dir)
        site_pdb = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb", dir=temp_dir)
        log_file = tempfile.NamedTemporaryFile(delete=False, dir=temp_dir)

        receptor.io.to_pdb(receptor_pdb.name, PDBOPTIONS)
        site.io.to_pdb(site_pdb.name, PDBOPTIONS)

        ligand_pdbs = []
        for lig in ligands:
            ligand_pdb = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb", dir=temp_dir)
            ComplexUtils.align_to(lig, receptor)
            lig.io.to_pdb(ligand_pdb.name, PDBOPTIONS)
            ligand_pdbs.append(ligand_pdb)

        smina_output_sdfs = []
        ligand_count = len(ligands)
        for ligand_pdb in ligand_pdbs:
            output_sdf = tempfile.NamedTemporaryFile(delete=False, prefix="output", suffix=".sdf", dir=temp_dir)
            process = self.run_smina(ligand_pdb, receptor_pdb, site_pdb, output_sdf, log_file, exhaustiveness, modes, autobox, ligand_count)
            self.handle_loading_bar(process, len(ligands))
            smina_output_sdfs.append(output_sdf)
        return smina_output_sdfs

    def run_smina(self, ligand_pdb, receptor_pdb, site_pdb, output_sdf, log_file, exhaustiveness=None, modes=None, autobox=None, ligand_count=1):
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

        Logs.debug("Run SMINA")
        cmd = [SMINA_PATH, *smina_args]
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        self.handle_loading_bar(process, ligand_count)
        return process

    def handle_loading_bar(self, process, ligand_count):
        """Render loading bar from stdout on the menu.

        stdout has a loading bar of asterisks. Every asterisk represents about 2% completed
        """
        stars_per_complex = 51
        total_stars = stars_per_complex * ligand_count

        for c in iter(lambda: process.stdout.read(1), b''):
            if c.decode() == '*':
                self.loading_bar_counter += 1
                self.plugin.update_loading_bar(self.loading_bar_counter, total_stars)
            sys.stdout.buffer.write(c)
