from Bio.PDB import PDBParser, Select
import numpy as np

def find_centroid_and_max_lengths(pdb_file):
    """
    Finds the centroid and maximum lengths along the xyz axes for a PDB structure.

    Args:
        pdb_file (str): The path to the PDB file.

    Returns:
        tuple: A tuple containing the centroid (x, y, z) and maximum lengths (x, y, z).
    """

    parser = PDBParser()
    structure = parser.get_structure("structure", pdb_file)

    class AtomSelector(Select):
        def accept_atom(self, atom):
            # You can customize the selection based on your needs (e.g., only atoms from a specific chain, residue type)
            return True

    selector = AtomSelector()

    # Convert the generator to a list
    selected_atoms = list(structure.get_atoms())

    # Calculate centroid
    centroid = sum(np.array(atom.get_coord()) for atom in selected_atoms) / len(selected_atoms)

    # Calculate maximum and minimum coordinates along each axis
    max_x = max(atom.get_coord()[0] for atom in selected_atoms)
    min_x = min(atom.get_coord()[0] for atom in selected_atoms)
    max_y = max(atom.get_coord()[1] for atom in selected_atoms)
    min_y = min(atom.get_coord()[1] for atom in selected_atoms)
    max_z = max(atom.get_coord()[2] for atom in selected_atoms)
    min_z = min(atom.get_coord()[2] for atom in selected_atoms)

    max_length_x = max_x - min_x
    max_length_y = max_y - min_y
    max_length_z = max_z - min_z

    return centroid, (max_length_x, max_length_y, max_length_z)

# Example usage:
pdb_file = "base_proteasome.pdb"
centroid, max_lengths = find_centroid_and_max_lengths(pdb_file)

print("Centroid:", centroid)
print("Maximum lengths along xyz axes:")
print("X:", max_lengths[0])
print("Y:", max_lengths[1])
print("Z:", max_lengths[2])
