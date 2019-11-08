import nanome

class ComplexUtils:

    @staticmethod
    def convert_atoms_to_absolute_position(complex):
        mat = complex.get_complex_to_workspace_matrix()
        for atom in complex.atoms:
            atom.position = mat * atom.position

    @staticmethod
    def convert_atoms_to_relative_position(complex, mat):
        # mat = reference.get_workspace_to_complex_matrix()
        for atom in complex.atoms:
            atom.position = mat * atom.position

    @staticmethod
    def combine_ligands(receptor, ligands):
        combined_ligands = nanome.structure.Complex()
        combined_ligands.names = []
        for ligand in ligands:
                ComplexUtils.convert_atoms_to_absolute_position(ligand)
                ComplexUtils.convert_atoms_to_relative_position(ligand, receptor.m_workspace_to_complex)
                combined_ligands.names.append(ligand.name)
                for molecule in ligand.molecules:
                    combined_ligands.add_molecule(molecule)
        return combined_ligands