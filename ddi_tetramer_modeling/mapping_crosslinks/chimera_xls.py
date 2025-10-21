# Load necessary ChimeraX modules
from chimerax.core.commands import run
from chimerax.core.errors import UserError
from chimerax.atomic import Residue, Atoms, AtomicStructure
from chimerax.geometry import distance

# Define the file paths (hardcoded)
PDB_FILE = "./input_data/Q8WTU0_V1_3.pdb"
#CROSSLINK_FILE = "pymol_data/formatted_random_lysine_doubles.csv"
CROSSLINK_FILE = "pymol_data/bifunc_pairs.csv"

# Function to load PDB structure
def load_structure(session, pdb_file):
    run(session, f'open {pdb_file}')

# Function to read crosslink data
def read_crosslink_data(file_path):
    crosslinks = []
    with open(file_path, 'r') as file:
        next(file)  # Skip the header line
        for line in file:
            if line.strip():
                chain1, res1, atom1, chain2, res2, atom2 = line.strip().split(',')
                crosslinks.append((chain1.strip(), int(res1.strip()), atom1.strip(), 
                                   chain2.strip(), int(res2.strip()), atom2.strip()))
    return crosslinks

# Function to get an atom by chain, residue number, and atom name
def find_atom_by_chain_residue(model, chain, residue_num, atom_name):
    for residue in model.residues:
        if residue.chain_id == chain and residue.number == residue_num:
            for atom in residue.atoms:
                if atom.name == atom_name:
                    return atom
    return None

# Function to visualize crosslinks on the structure
def plot_crosslinks(session, pdb_model, crosslinks):
    violation_count = 0
    processed_crosslinks = set()
    
    for crosslink in crosslinks:
        chain1, res1, atom1, chain2, res2, atom2 = crosslink

        # Skip duplicate crosslinks
        if (chain1, res1, atom1, chain2, res2, atom2) in processed_crosslinks or \
           (chain2, res2, atom2, chain1, res1, atom1) in processed_crosslinks:
            continue
        
        processed_crosslinks.add((chain1, res1, atom1, chain2, res2, atom2))

        # Find atoms in the structure by chain, residue number, and atom name
        atom1_obj = find_atom_by_chain_residue(pdb_model, chain1, res1, atom1)
        atom2_obj = find_atom_by_chain_residue(pdb_model, chain2, res2, atom2)
        
        # If both atoms are found, create or modify a link (distance line)
        if atom1_obj and atom2_obj:
            dist = distance(atom1_obj.coord, atom2_obj.coord)
            color = 'dark green' if dist < 30 else 'dark red'
            if dist >= 30:
                violation_count += 1
            try:
                run(session, f'distance {atom1_obj.atomspec} {atom2_obj.atomspec} color {color} radius 0.6')
            except UserError:
                run(session, f'distance style {atom1_obj.atomspec} {atom2_obj.atomspec} color {color} radius 0.6')
    
    print(f'Number of crosslink distance violations: {violation_count}')

# Main function to run in ChimeraX
def visualize_crosslinks(session):
    # Load the PDB structure
    load_structure(session, PDB_FILE)
    
    # Get the currently opened atomic structure model
    pdb_model = None
    for model in session.models:
        if isinstance(model, AtomicStructure):
            pdb_model = model
            break
    
    if pdb_model is None:
        raise ValueError("No atomic structure model found.")
    
    # Apply cartoon representation and color to the protein
    run(session, 'cartoon')
    run(session, 'color khaki cartoons')
    
    # Read the crosslink data
    crosslinks = read_crosslink_data(CROSSLINK_FILE)
    
    # Plot crosslinks on the structure
    plot_crosslinks(session, pdb_model, crosslinks)

# Example of usage:
# In ChimeraX: visualize_crosslinks(session)
