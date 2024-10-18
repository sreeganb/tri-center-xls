#!/usr/bin/env python3.11

import pandas as pd
import os

# Define the mapping from Protein to Chain (corrected)
protein_to_chain = {
    'Rpt6': 'x',
    'Rpt3': 'y',
    'Rpt4': 'z',
    'Rpn2': '1',
    'Rpt5': '0',
    'Rpt2': 'w',
    'Rpt1': 'v'
}

def ensure_directory(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

def format_for_pymol(input_file, output_file):
    # Check if the input file exists
    if not os.path.exists(input_file):
        print(f"Error: The file {input_file} does not exist.")
        return

    # Read the input file into a DataFrame
    df = pd.read_csv(input_file)

    # Use the mapping in protein_to_chain to add columns next to Protein1 and Protein2 as chain1 and chain2
    df['chain1'] = df['Protein1'].map(protein_to_chain)
    df['chain2'] = df['Protein2'].map(protein_to_chain)
    
    # Add two columns 'Atom1' and 'Atom2' with value 'CA'
    df['Atom1'] = 'CA'
    df['Atom2'] = 'CA'
    
    # Drop the columns 'Protein1' and 'Protein2' and rearrange the columns in the order: chain1, Residue1, Atom1, chain2, Residue2, Atom2
    df = df[['chain1', 'Residue1', 'Atom1', 'chain2', 'Residue2', 'Atom2']]
    
    # Print the first few rows of the DataFrame
    print(df.head())
    
    # Write the DataFrame to the output file
    df.to_csv(output_file, index=False)

def main():
    # Ensure the output directory exists
    output_directory = 'pymol_data'
    ensure_directory(output_directory)
    
    # Format random lysine doubles for PyMOL
    format_for_pymol('filtered_cleaned_double_links.csv', os.path.join(output_directory, 'dsso-xls.csv'))

if __name__ == "__main__":
    main()