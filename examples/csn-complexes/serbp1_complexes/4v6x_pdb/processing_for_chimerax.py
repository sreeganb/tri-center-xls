#!/usr/bin/env python3.11

import pandas as pd
import os

# Define the mapping from Protein to Chain (corrected)
protein_to_chain = {
    'RPS12':  'G',
    'SERBP1': 'S',
    'RPL27':  'K',
    'RPL34':  'I',
    'RPL7A':  'J',
    'RPS10':  'E',
    'RPS11':  'F',
    'RPS12':  'G',
    'RPS30':  'C'
}

def ensure_directory(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

def process_triple_links(input_file, output_file):
    # Read the CSV file into a DataFrame, skipping the first row
    headers = ['Protein1', 'Residue1', 'Protein2', 'Residue2', 'Protein3', 'Residue3']
    df = pd.read_csv(input_file, header=0, names=headers)

    # Initialize an empty list to store the new rows
    new_rows = []

    # Process each row in the DataFrame
    for index, row in df.iterrows():
        new_rows.append([row['Protein1'], row['Residue1'], row['Protein2'], row['Residue2']])
        new_rows.append([row['Protein2'], row['Residue2'], row['Protein3'], row['Residue3']])
        new_rows.append([row['Protein3'], row['Residue3'], row['Protein1'], row['Residue1']])

    # Create a new DataFrame from the new rows
    new_df = pd.DataFrame(new_rows, columns=['Protein1', 'Residue1', 'Protein2', 'Residue2'])
    new_df.to_csv(output_file, index=False)

def process_regular_xls(input_file, output_file):
    # Read the CSV file into a DataFrame, skipping the first row
    headers = ['Protein1', 'Residue1', 'Protein2', 'Residue2'] 
    df = pd.read_csv(input_file, header=0, names=headers)

    # No need for further processing, just save to output file
    df.to_csv(output_file, index=False)


def format_for_pymol(input_file, output_file):
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
    output_directory = 'data_chimerax'
    ensure_directory(output_directory)

    # Set this variable to 1 to process double crosslinks, 0 for triple crosslinks
    doubles = 0  

    if doubles:
        process_regular_xls('serbp1_bifunctional.csv', os.path.join(output_directory, 'paired_double_links_serbp1.csv'))
        format_for_pymol(os.path.join(output_directory, 'paired_double_links_serbp1.csv'), os.path.join(output_directory, 'formatted_doubles_serbp1.csv'))
    else:
        # Process triple links
        process_triple_links('serbp1_trifunctional.csv', os.path.join(output_directory, 'paired_triple_links_serbp1.csv'))
        format_for_pymol(os.path.join(output_directory, 'paired_triple_links_serbp1.csv'), os.path.join(output_directory, 'formatted_triplets_serbp1.csv'))

if __name__ == "__main__":
    main()
