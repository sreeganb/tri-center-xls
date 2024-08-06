#!/usr/bin/env python3.11

import pandas as pd
import re
from Bio import PDB, SeqIO
import warnings
import numpy as np
import csv
import os
import matplotlib.pyplot as plt
import seaborn as sns
from Bio.PDB import MMCIFParser, PDBIO, Select


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

def replace_alpha_beta(input_file, output_file):
    """Replaces α and β with alpha and beta in specified columns."""
    alpha_beta_pattern = r"α(\d+)|β(\d+)"

    with open(input_file, mode='r') as infile, open(output_file, mode='w', newline='') as outfile:
        reader = csv.DictReader(infile)
        fieldnames = reader.fieldnames
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()

        for row in reader:  

            for column in ['Common/Gene Name A', 'Common/Gene Name B', 'Common/Gene Name C']:
                match = re.search(alpha_beta_pattern, row[column])
                if match:
                    if match.group(1):
                        row[column] = f"alpha{match.group(1)}"
                    else:
                        row[column] = f"beta{match.group(2)}"
            writer.writerow(row)

def extract_protein_info(fasta_file, csv_file):
    """Extracts protein name, chain ID, and common name from FASTA and CSV files into a DataFrame."""
    data = []
    auth_pattern = r"\[auth (\w+)\]"

    # Read CSV file
    csv_df = pd.read_csv(csv_file)

    with open(fasta_file, "r") as f_in, open("mod_5gjr.fasta", "w") as f_out:
        for record in SeqIO.parse(f_in, "fasta"):
            header = record.description
            parts = header.split("|")

            protein_name = parts[2] if len(parts) > 2 else "Unknown"
            chain_id = ','.join(re.findall(auth_pattern, header)) if re.findall(auth_pattern, header) else "Unknown"

            # Find corresponding common name from CSV
            common_name = csv_df[csv_df['Subunit Name A'] == protein_name]['Common/Gene Name A'].values
            common_name = common_name[0] if len(common_name) > 0 else "Unknown"

            data.append({'protein_name': protein_name, 'chain_id': chain_id, 'common_name': common_name})

            # Write modified FASTA record to new file
            f_out.write(f">{common_name}\n")
            f_out.write(str(record.seq) + "\n")

    return pd.DataFrame(data)

def extract_base_of_proteasome(input_file, output_file):
    """Extracts rows where 'Common/Gene Name A' and 'Common/Gene Name B' are from base subunits,
    and 'Common/Gene Name C' is either empty or from base subunits, and writes to a new CSV file."""
    # Define the list of subunit names to filter
    base_subunits = {'Rpt1', 'Rpt2', 'Rpt3', 'Rpt4', 'Rpt5', 'Rpt6', 'Rpn1', 'Rpn2', 'Rpn10', 'Rpn13'}

    # Read the input CSV file
    df = pd.read_csv(input_file)

    # Check if 'Common/Gene Name A' and 'Common/Gene Name B' are in base_subunits
    # and 'Common/Gene Name C' is either empty or in base_subunits
    mask = df.apply(
        lambda row: all(
            value in base_subunits for value in [row['Common/Gene Name A'], row['Common/Gene Name B']]
        ) and (pd.isna(row['Common/Gene Name C']) or row['Common/Gene Name C'] in base_subunits),
        axis=1
    )

    # Filter rows based on the mask
    filtered_df = df[mask]

    # Write the filtered DataFrame to a new CSV file
    filtered_df.to_csv(output_file, index=False)

def extract_all_chains_to_pdb(mmcif_file, output_pdb):
    """Extracts all chains from an mmCIF file and combines them into a single PDB file.

    Args:
        mmcif_file (str): Path to the input mmCIF file.
        output_pdb (str): Path to the output PDB file.
    """

    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("structure", mmcif_file)

    io = PDBIO()
    counter = 0
    for model in structure:
        for chain in model:
            chain_id = chr(ord('A') + counter)  # Assign a new chain ID
            counter += 1
            io.set_structure(chain)
            io.save(output_pdb, select=lambda x: x.id == chain_id, append=True)

def main():
    fasta_file = "../fasta/rcsb_pdb_5GJR.fasta"
    fasta_dir = "../fasta/"
    pdb_dir_1 = "./"
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
    
    remove_duplicate_rows('processed-data.csv', 'processed-data_1.csv')
    replace_alpha_beta('processed-data_1.csv', 'unique-xls.csv')
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
    # Call the function with the appropriate file paths
    input_file = 'new-unique-xls.csv'
    output_file = 'base-of-proteasome.csv'
    extract_base_of_proteasome(input_file, output_file)
    seqdat = extract_protein_info(fasta_dir + 'rcsb_pdb_5GJR.fasta', pdb_dir_1 + 'new-unique-xls.csv')
    
    df1 = df1_exploded
    
    # Using the newly created dataframe for the base of the proteasome here
    # We should check with this and then revert back to the full proteasome 
    # or not, doesnt matter for now at least
    base_df = pd.read_csv('./base-of-proteasome.csv')
    df1 = base_df
    # Plot the density of the proteins that are crosslinked
    comb_data = pd.concat([df1['Common/Gene Name A'], df1['Common/Gene Name B'], df1['Common/Gene Name C']])
    value_counts = comb_data.value_counts().sort_index()

    # Create a density plot with seaborn
    sns.set_style("whitegrid")  # Set seaborn style for better aesthetics
    plt.figure(figsize=(10, 7))
    sns.histplot(comb_data, bins=30, kde=False)  # Use histplot for histogram
    #sns.kdeplot(value_counts, shade=True, color="skyblue")  # Use kernel density estimation
    
    # Customize the plot
    plt.xlabel('Unique Crosslinked Proteins')
    plt.ylabel('Frequency')
    plt.title('Frequency of Crosslinked Proteins')
    plt.xticks(rotation=45)  # Rotate x-axis labels for readability
    plt.savefig('xl_density.pdf')
    plt.tight_layout()
    plt.show()

    cif_file = './5gjr.cif'
    # Create a list of unique elements
    # Split the entries and create separate lists for first and second elements
    first_elements = []
    second_elements = []

    for entry in df1['chain ID A']:
        first, second = entry.split(',')
        first_elements.append(first)
        second_elements.append(second)

    # Remove duplicates by converting to sets and then back to lists
    unique_first_elements = list(set(first_elements))
    unique_second_elements = list(set(second_elements))

    # Print the results
    print("Unique first elements:", unique_first_elements)
    print("Unique second elements:", unique_second_elements)

    #unique_chain_ids = list(set(sum(df1['chain ID A'].str.split(','), [])))
    #print("the unique chain IDs are: ", unique_chain_ids)
    parser = PDB.MMCIFParser()
    structure = parser.get_structure('structure', cif_file)
    
    calpha_distances = []
    calpha_doubles = []
    non_empty_entries = df1['XL C'].notnull().sum()
    print(f"Number of non-empty entries: {non_empty_entries}")
    
    not_found = 0
    valid_entries = []  # To store valid residue pairs and chain IDs
    n_double_links = len(df1['XL A'])  # Number of double links
    double_links = []
    triple_links = []
    for i in range(n_double_links):
        if not df1.loc[i, 'XL A'].startswith('K') or not df1.loc[i, 'XL B'].startswith('K'):
            continue
        res1_num = int(df1.loc[i, 'XL A'][1:])  # Drop the first character and convert to integer
        res2_num = int(df1.loc[i, 'XL B'][1:])  # Drop the first character and convert to integer
        
        chain1_id = df1.at[i, 'chain ID A'].split(',')[0]
        chain2_id = df1.at[i, 'chain ID B'].split(',')[0]
        dis = calculate_ca_distance(structure, res1_num, chain1_id, res2_num, chain2_id)
        if dis != -1:
            dis_alt = 99999
            if dis > 35.0:
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
            #calpha_distances.append(dis)
        else:
            not_found += 1
        
        if min(dis,dis_alt) < 35.0:
            valid_entries.append([res1_num, chain1_id, res2_num, chain2_id])
        if i > non_empty_entries-1: # Only process double links
            double_links.append([res1_num, chain1_id, res2_num, chain2_id])
            calpha_doubles.append(min(dis, dis_alt))
        if i < non_empty_entries: # Only process triple links
            calpha_distances.append(min(dis, dis_alt))
            triple_links.append([res1_num, chain1_id, res2_num, chain2_id])
                  
        if i < non_empty_entries:
            # Only process triple links
            if not df1.loc[i, 'XL A'].startswith('K') or not df1.loc[i, 'XL B'].startswith('K') or not df1.loc[i, 'XL C'].startswith('K'):
                continue
            res3_num = int(df1.loc[i, 'XL C'][1:])  # Drop the first character and convert to integer
            chain3_id = df1.at[i, 'chain ID C'].split(',')[0]
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
            if min(dis, dis_alt) < 35.0:
                valid_entries.append([res1_num, chain1_id, res3_num, chain3_id])
                triple_links.append([res1_num, chain1_id, res3_num, chain3_id])
                    
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
            if min(dis, dis_alt) < 35.0:
                valid_entries.append([res2_num, chain2_id, res3_num, chain3_id])
                triple_links.append([res2_num, chain2_id, res3_num, chain3_id])
                    
    print(f"Number of not found entries: {not_found}")
    # Save valid entries to a CSV file
    valid_df = pd.DataFrame(valid_entries, columns=['Residue 1 Number', 'Chain 1 ID', 'Residue 2 Number', 'Chain 2 ID'])
    valid_df.to_csv('valid_distances.csv', index=False)
    double_df = pd.DataFrame(double_links, columns=['Residue 1 Number', 'Chain 1 ID', 'Residue 2 Number', 'Chain 2 ID'])
    double_df.to_csv('double_links.csv', index=False)
    triple_df = pd.DataFrame(triple_links, columns=['Residue 1 Number', 'Chain 1 ID', 'Residue 2 Number', 'Chain 2 ID'])
    triple_df.to_csv('triple_links.csv', index=False)
    
    #print(f"Number of valid entries: {len(calpha_distances)}")
    print(f"Number of double links: {len(calpha_doubles)}")
    print(f"Number of triple links: {len(calpha_distances)}")

    plt.hist(calpha_distances, edgecolor='black', alpha=0.3, label='triple links')
    plt.hist(calpha_doubles, edgecolor='black', alpha=0.3, label='double links')
    plt.axvline(x=35.0, color='red', linestyle='--')
    plt.xlabel('C-alpha Distance (Å)')
    plt.ylabel('Frequency')
    plt.legend()
    plt.title(r"Distribution of C$\alpha$ distances of XL Lysine residues")
    plt.tight_layout()
    plt.savefig('triple_links.pdf')
    plt.show()
    
if __name__ == "__main__":
    main()