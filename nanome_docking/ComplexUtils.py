import nanome

class ComplexUtils:

    @staticmethod
    def align_to(complex, reference_complex):
        m = complex.get_complex_to_workspace_matrix()
        for atom in complex.atoms:
            atom.position = m * atom.position
        complex.position = reference_complex.position
        complex.rotation = reference_complex.rotation
        m = complex.get_workspace_to_complex_matrix()
        for atom in complex.atoms:
            atom.position = m * atom.position

    @staticmethod
    def combine_ligands(receptor, ligands):
        combined_ligands = nanome.structure.Complex()
        combined_ligands.names = []
        for ligand in ligands:
            ComplexUtils.align_to(ligand, receptor)
            combined_ligands.names.append(ligand.full_name)
            for molecule in ligand.molecules:
                combined_ligands.add_molecule(molecule)
        return combined_ligands
