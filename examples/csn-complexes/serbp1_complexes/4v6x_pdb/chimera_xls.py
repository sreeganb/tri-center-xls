import os
import logging

# Load necessary ChimeraX modules
from chimerax.core.commands import run
from chimerax.core.errors import UserError
from chimerax.atomic import Residue, Atoms, AtomicStructure
from chimerax.geometry import distance
import numpy as np

# Define the file paths (hardcoded)
PDB_FILE = "combined_chains_4v6x.pdb"
#CROSSLINK_FILE = "data_chimerax/formatted_triplets_serbp1.csv"
CROSSLINK_FILE = "data_chimerax/formatted_doubles_serbp1.csv"

# Configure logging
LOG_DIR = "chimerax_logs"
if not os.path.exists(LOG_DIR):
    os.makedirs(LOG_DIR)
LOG_FILE = os.path.join(LOG_DIR, "crosslink_analysis.log")

logging.basicConfig(filename=LOG_FILE, level=logging.INFO, 
                    format='%(asctime)s - %(levelname)s - %(message)s')

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
    distances = []
    missing_residues = set()
    
    # Get all residue numbers from the PDB model
    pdb_residues = {(residue.chain_id, residue.number) for residue in pdb_model.residues}

    for crosslink in crosslinks:
        chain1, res1, atom1, chain2, res2, atom2 = crosslink

        # Check for missing residues
        if (chain1, res1) not in pdb_residues:
            missing_residues.add((chain1, res1))
        if (chain2, res2) not in pdb_residues:
            missing_residues.add((chain2, res2))

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
            color = 'dark green' if dist < 32 else 'dark red'
            if dist >= 32:
                violation_count += 1
            try:
                run(session, f'distance {atom1_obj.atomspec} {atom2_obj.atomspec} color {color} radius 0.6')
            except UserError:
                run(session, f'distance style {atom1_obj.atomspec} {atom2_obj.atomspec} color {color} radius 0.6')
            distances.append(dist)

    # --- Output and Logging ---
    log_and_print(f"Number of crosslinks read from file: {len(crosslinks)}")
    log_and_print(f"Number of crosslink distance violations: {violation_count}")
    log_and_print(f"Mean distance: {np.mean(distances):.2f} Å")
    log_and_print(f"Standard deviation: {np.std(distances):.2f} Å")
    log_and_print(f"Maximum distance: {np.max(distances):.2f} Å")
    log_and_print(f"Minimum distance: {np.min(distances):.2f} Å")
    log_and_print(f"Number of residues in crosslink data not found in PDB: {len(missing_residues)}")
    if missing_residues:
        log_and_print(f"Missing residues: {missing_residues}")

    # Additional structural information
    log_and_print(f"Number of chains in PDB: {len(pdb_model.chains)}")
    log_and_print(f"Total number of residues in PDB: {len(pdb_model.residues)}")

# Helper function to log and print
def log_and_print(message):
    print(message)
    logging.info(message)

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
