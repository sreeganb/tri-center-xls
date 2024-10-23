#!/usr/bin/env python3.11

import pandas as pd
import os
from Bio import PDB
import numpy as np
import matplotlib.pyplot as plt

# Define the mapping from Protein to Chain (corrected)
# trying to be exhaustive here, get all the names possibly 
# from the experimental dataset
protein_to_chain = {
    'Rpt6': 'x',
    'Rpt3': 'y',
    'Rpt4': 'z',
    'Rpn2': '1',
    'Rpt5': '0',
    'Rpt2': 'w',
    'Rpt1': 'v',
    'Rpn11': '9',
    'Rpn3': '6',
    'Rpn5': '3',
    'Rpn9': '2',
    'alpha3': 'D',
    'alpha6': 'G',
    'alpha7': 'X',
    'beta5': 'F',
    'alpha4': 'E',
    'alpha1': 'B',
    'beta6': 'f',
    'alpha2': 'C',
    'beta4': 'd'
}

def ensure_directory(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

def get_atom_coordinates(pdb_file, chain, residue, atom):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('structure', pdb_file)
    for model in structure:
        for chain_obj in model:
            if chain_obj.id == chain:
                for residue_obj in chain_obj:
                    if residue_obj.id[1] == residue:
                        for atom_obj in residue_obj:
                            if atom_obj.id == atom:
                                return atom_obj.coord
    return None

def calculate_distance(coord1, coord2):
    return np.linalg.norm(coord1 - coord2)

def format_data(input_file, output_file, pdb_file):
    # Read the input file into a DataFrame
    # data in the form of Protein1, Residue1, Protein2, Residue2
    df = pd.read_csv(input_file)

    # Use the mapping in protein_to_chain to add columns next to Protein1 and Protein2 as chain1 and chain2
    df['chain1'] = df['Protein1'].map(protein_to_chain)
    df['chain2'] = df['Protein2'].map(protein_to_chain)
    
    # Add two columns 'Atom1' and 'Atom2' with value 'CA'
    df['Atom1'] = 'CA'
    df['Atom2'] = 'CA'
    
    # Calculate distances between Atom1 and Atom2
    distances = []
    for index, row in df.iterrows():
        coord1 = get_atom_coordinates(pdb_file, row['chain1'], row['Residue1'], row['Atom1'])
        coord2 = get_atom_coordinates(pdb_file, row['chain2'], row['Residue2'], row['Atom2'])
        if coord1 is not None and coord2 is not None:
            distance = calculate_distance(coord1, coord2)
            distances.append(distance)
        else:
            distances.append(None)
    
    df['Distance'] = distances
    
    # Drop the columns 'Protein1' and 'Protein2' and rearrange the columns in the order: chain1, Residue1, Atom1, chain2, Residue2, Atom2, Distance
    df = df[['chain1', 'Residue1', 'Atom1', 'chain2', 'Residue2', 'Atom2', 'Distance']]
    
    # Write the DataFrame to the output file
    df.to_csv(output_file, index=False)
    
    # Plot the distance distribution
    plot_distance_distribution(distances, 'output_data')

def plot_distance_distribution(distances, output_directory):
    # Remove None values from distances
    distances = [d for d in distances if d is not None]
    
    plt.hist(distances, bins=30, edgecolor='black')
    plt.title('Distribution of Distances Between C-alpha Atoms of Lysine Residues')
    plt.xlabel('Distance (Å)')
    plt.ylabel('Frequency')
    
    # Save the plot to the output directory
    plt.savefig(os.path.join(output_directory, 'distance_distribution.png'))
    plt.show()

def main():
    # Ensure the output directory exists
    output_directory = 'output_data'
    ensure_directory(output_directory)
    
    # Format random lysine doubles and triplets for PyMOL
    format_data('exp_triples.csv', os.path.join(output_directory, 'exp_triples_distances.csv'), 'base_proteasome.pdb')
    
if __name__ == "__main__":
    main()