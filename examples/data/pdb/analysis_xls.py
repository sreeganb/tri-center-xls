#!/usr/bin/env python3.11

import pandas as pd
import re
from Bio import PDB, SeqIO
import warnings
import numpy as np
import csv
import os
import matplotlib.pyplot as plt

# Suppress all warnings
warnings.filterwarnings("ignore")

def extract_protein_info_to_dataframe(fasta_file):
    """Extracts protein name and chain ID from a FASTA file into a DataFrame."""
    data = []
    auth_pattern = r"\[auth (\w+)\]"
    
    for record in SeqIO.parse(fasta_file, "fasta"):
        header = record.description
        parts = header.split("|")
        
        protein_name = parts[2] if len(parts) > 2 else "Unknown"
        chain_id = ','.join(re.findall(auth_pattern, header)) if re.findall(auth_pattern, header) else "Unknown"
        
        data.append({'protein_name': protein_name, 'chain_id': chain_id})
    
    return pd.DataFrame(data)

def match_and_add_chain_ids(df, psome_df):
    """Matches proteins in df with those in psome_df and adds chain ID columns."""
    for col in ['chain ID A', 'chain ID B', 'chain ID C']:
        psome_df[col] = ''

    for index, row in psome_df.iterrows():
        protein_a = str(row['Subunit Name A']).lower()
        protein_b = str(row['Subunit Name B']).lower()
        protein_c = str(row['Subunit Name C']).lower()

        for _, df_row in df.iterrows():
            if df_row['protein_name'].lower() == protein_a:
                psome_df.at[index, 'chain ID A'] = df_row['chain_id']
            if df_row['protein_name'].lower() == protein_b:
                psome_df.at[index, 'chain ID B'] = df_row['chain_id']
            if df_row['protein_name'].lower() == protein_c:
                psome_df.at[index, 'chain ID C'] = df_row['chain_id']
    
    return psome_df

def remove_rows_with_empty_chain_id_a_b(df):
    """Removes rows where 'chain ID A' or 'chain ID B' is empty."""
    return df[(df['chain ID A'] != '') & (df['chain ID B'] != '')]

def remove_duplicate_rows(input_file, output_file):
    """Removes duplicate rows from the CSV file based on Subunit columns."""
    unique_rows = {}
    
    with open(input_file, mode='r') as infile:
        reader = csv.DictReader(infile)
        fieldnames = reader.fieldnames
        
        for row in reader:
            key = (row['Subunit A'], row['Subunit B'], row['Subunit C'])
            if key not in unique_rows:
                unique_rows[key] = row
    
    with open(output_file, mode='w', newline='') as outfile:
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()
        for row in unique_rows.values():
            writer.writerow(row)

def calculate_ca_distance(structure, res1_num, chain1_id, res2_num, chain2_id):
    """Calculates the distance between two C-alpha atoms in a structure."""
    ca_res1 = None
    ca_res2 = None

    for model in structure:
        for chain in model:
            if chain.id == chain1_id:
                for residue in chain:
                    if residue.get_id()[1] == res1_num and 'CA' in residue:
                        ca_res1 = residue['CA']
                        if residue.get_resname() != 'LYS':
                            print(f"The specified residue 1 is not LYS but instead it is {residue.get_resname()}")
                            return -1
            if chain.id == chain2_id:
                for residue in chain:
                    if residue.get_id()[1] == res2_num and 'CA' in residue:
                        ca_res2 = residue['CA']
                        if residue.get_resname() != 'LYS':
                            print(f"The specified residue 2 is not LYS but instead it is {residue.get_resname()}")
                            return -1
                          
    if ca_res1 is None or ca_res2 is None:
        print("Error: One or both specified C-alpha atoms were not found.")
        return -1

    return np.linalg.norm(ca_res1.get_coord() - ca_res2.get_coord())

def main():
    fasta_file = "../fasta/rcsb_pdb_5GJR.fasta"
    protein_df = extract_protein_info_to_dataframe(fasta_file)
    protein_df['chain_id'] = protein_df['chain_id'].astype(str)
    protein_df.to_csv('output.csv', index=False)
    
    psome_file = '../../derived_data/xl/proteasome-tri-crosslinks.csv'
    psome_df = pd.read_csv(psome_file)
    print(psome_df.head())
    
    # Initialize the 'chain ID' columns
    for col in ['chain ID A', 'chain ID B', 'chain ID C']:
        psome_df[col] = ''
    
    processed_data = psome_df.iloc[:, :12]  # Select the first 12 columns
    processed_data['Average Score'] = psome_df[['Score A', 'Score B', 'Score C']].mean(axis=1)
    
    processed_data = match_and_add_chain_ids(protein_df, processed_data)
    
    processed_data = remove_rows_with_empty_chain_id_a_b(processed_data)
    
    processed_data.to_csv('processed-data.csv', index=False)
    
    remove_duplicate_rows('processed-data.csv', 'unique-xls.csv')
    df1 = pd.read_csv('unique-xls.csv')
    
    # Functionality to split data in 'XL A', 'XL B', and 'XL C' columns   
    df1['XL A'] = df1['XL A'].str.split(r'[;|]')
    df1['XL B'] = df1['XL B'].str.split(r'[;|]')
    df1['XL C'] = df1['XL C'].str.split(r'[;|]')
    
    # Explode the 'Values' column to separate rows
    df1_exploded = df1.explode('XL A').explode('XL B').explode('XL C')
    # Reset index if needed
    df1_exploded.reset_index(drop=True, inplace=True)
    df1_exploded.to_csv('new-unique-xls.csv', index=False)
    
    df1 = df1_exploded
    cif_file = './5gjr.cif'
    parser = PDB.MMCIFParser()
    structure = parser.get_structure('structure', cif_file)
    
    calpha_distances = []
    non_empty_entries = df1['XL C'].notnull().sum()
    print(f"Number of non-empty entries: {non_empty_entries}")
    
    not_found = 0
    valid_entries = []  # To store valid residue pairs and chain IDs
    n_double_links = len(df1['XL A'])  # Number of double links
    print(f"Number of double links: {n_double_links}")
    for i in range(non_empty_entries):
        if not df1.loc[i, 'XL A'].startswith('K') or not df1.loc[i, 'XL B'].startswith('K') or not df1.loc[i, 'XL C'].startswith('K'):
            continue
        res1_num = int(df1.loc[i, 'XL A'][1:])  # Drop the first character and convert to integer
        res2_num = int(df1.loc[i, 'XL B'][1:])  # Drop the first character and convert to integer
        res3_num = int(df1.loc[i, 'XL C'][1:])  # Drop the first character and convert to integer
            
        chain1_id = df1.at[i, 'chain ID A'].split(',')[0]
        chain2_id = df1.at[i, 'chain ID B'].split(',')[0]
        chain3_id = df1.at[i, 'chain ID C'].split(',')[0]
            
        dis = calculate_ca_distance(structure, res1_num, chain1_id, res2_num, chain2_id)
        if dis != -1:
            dis_alt = 99999
            if dis > 35.0:
                #print(f"Distance between C-alpha atoms: {dis:.2f} Å and the chain ID is {chain1_id}, {chain2_id}")
                # Try with the second element if distance is high
                try:
                    chain1_id_alt = df1.at[i, 'chain ID A'].split(',')[1]
                    dis_alt = calculate_ca_distance(structure, res1_num, chain1_id_alt, res2_num, chain2_id)
                    if dis_alt != -1 and dis_alt < 35.0:
                        chain1_id = chain1_id_alt
                    #if dis_alt != -1:
                        #calpha_distances.append(dis_alt)
                        #print(f"Trying alternative chain ID: Distance: {dis_alt:.2f} Å, Chain IDs: {chain1_id_alt}, {chain2_id}")
                except IndexError:
                    # Handle cases where there's no second element
                    pass
            calpha_distances.append(min(dis, dis_alt))
            #calpha_distances.append(dis)
        else:
            not_found += 1
        if dis < 35.0:
            valid_entries.append([res1_num, chain1_id, res2_num, chain2_id])
        
        dis = calculate_ca_distance(structure, res1_num, chain1_id, res3_num, chain3_id)
        if dis != -1:
            dis_alt = 99999
            if dis > 35.0:
                print(f"Distance between C-alpha atoms: {dis:.2f} Å and the chain ID is {chain1_id}, {chain2_id}")
                # Try with the second element if distance is high
                try:
                    chain1_id_alt = df1.at[i, 'chain ID A'].split(',')[1]
                    dis_alt = calculate_ca_distance(structure, res1_num, chain1_id_alt, res3_num, chain3_id)
                    if dis_alt != -1 and dis_alt < 35.0:
                        chain1_id = chain1_id_alt
                        #calpha_distances.append(dis_alt)
                        print(f"Trying alternative chain ID: Distance: {dis_alt:.2f} Å, Chain IDs: {chain1_id_alt}, {chain2_id}")
                except IndexError:
                    # Handle cases where there's no second element
                    pass
            calpha_distances.append(min(dis, dis_alt))
        else:
            not_found += 1
        if dis < 35.0:
            valid_entries.append([res1_num, chain1_id, res3_num, chain3_id])
                    
        dis = calculate_ca_distance(structure, res2_num, chain2_id, res3_num, chain3_id)
        if dis != -1:
            dis_alt = 99999
            if dis > 35.0:
                print(f"Distance between C-alpha atoms: {dis:.2f} Å and the chain ID is {chain2_id}, {chain3_id}")
                # Try with the second element if distance is high
                try:
                    chain2_id_alt = df1.at[i, 'chain ID B'].split(',')[1]
                    dis_alt = calculate_ca_distance(structure, res1_num, chain2_id_alt, res3_num, chain3_id)
                    if dis_alt != -1 and dis_alt < 35.0:
                        #calpha_distances.append(dis_alt)
                        chain2_id = chain2_id_alt
                        print(f"Trying alternative chain ID: Distance: {dis_alt:.2f} Å, Chain IDs: {chain2_id_alt}, {chain3_id}")
                except IndexError:
                    # Handle cases where there's no second element
                    pass
            calpha_distances.append(min(dis, dis_alt))
            #calpha_distances.append(dis)
        else:
            not_found += 1
        if dis < 35.0:
            valid_entries.append([res2_num, chain2_id, res3_num, chain3_id])
        
    print(f"Number of not found entries: {not_found}")
    # Save valid entries to a CSV file
    valid_df = pd.DataFrame(valid_entries, columns=['Residue 1 Number', 'Chain 1 ID', 'Residue 2 Number', 'Chain 2 ID'])
    valid_df.to_csv('valid_distances.csv', index=False)

    plt.hist(calpha_distances, bins=20, edgecolor='black')
    plt.axvline(x=35.0, color='red', linestyle='--')
    plt.xlabel('C-alpha Distance (Å)')
    plt.ylabel('Frequency')
    plt.title(r"Distribution of C$\alpha$ distances of trilinks")
    plt.tight_layout()
    plt.savefig('triple_links.pdf')
    plt.show()
    
if __name__ == "__main__":
    main()
