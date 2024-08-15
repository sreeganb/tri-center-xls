#!/usr/bin/env python
import pandas as pd
from Bio.PDB import PDBParser
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import re

def filter_chain_ids(df, alphabets):
    """
    Filters the chain ID columns based on the provided list of alphabets.
    
    Args:
        df (pd.DataFrame): The input dataframe containing chain ID columns.
        alphabets (list[str]): The list of alphabets to check against.

    Returns:
        pd.DataFrame: The modified dataframe with filtered chain IDs.
    """
    chain_columns = ['chain ID A', 'chain ID B', 'chain ID C']
    
    def filter_entry(entry):
        if pd.isna(entry):
            return entry
        parts = entry.split(',')
        filtered_parts = [part for part in parts if part in alphabets]
        return ','.join(filtered_parts) if filtered_parts else ''

    for column in chain_columns:
        df[column] = df[column].apply(filter_entry)
    
    return df

def calculate_ca_distance(structure, res1_num, chain1_id, res2_num, chain2_id):
    """Calculates the distance between two C-alpha atoms in a structure."""
    ca_res1 = None
    ca_res2 = None

    # Initialize variables for missing information
    missing_residues = []
    missing_chains = []

    for model in structure:
        for chain in model:
            if chain.id == chain1_id:
                found_residue = False
                for residue in chain:
                    if residue.get_id()[1] == res1_num and 'CA' in residue:
                        ca_res1 = residue['CA']
                        if residue.get_resname() != 'LYS':
                            print(f"The specified residue 1 is not LYS but instead it is {residue.get_resname()}")
                            return np.nan
                        found_residue = True
                if not found_residue:
                    missing_residues.append((res1_num, chain1_id))
            if chain.id == chain2_id:
                found_residue = False
                for residue in chain:
                    if residue.get_id()[1] == res2_num and 'CA' in residue:
                        ca_res2 = residue['CA']
                        if residue.get_resname() != 'LYS':
                            print(f"The specified residue 2 is not LYS but instead it is {residue.get_resname()}")
                            return np.nan
                        found_residue = True
                if not found_residue:
                    missing_residues.append((res2_num, chain2_id))
                    
            if chain.id not in [chain1_id, chain2_id]:
                missing_chains.append(chain.id)

    if ca_res1 is None or ca_res2 is None:
        if ca_res1 is None:
            print(f"Error: C-alpha atom not found for residue {res1_num} in chain {chain1_id}.")
        if ca_res2 is None:
            print(f"Error: C-alpha atom not found for residue {res2_num} in chain {chain2_id}.")
        if missing_residues:
            print(f"Missing residues and chains: {missing_residues}")
        if missing_chains:
            print(f"Chains not present in structure: {missing_chains}")
        return np.nan

    return np.linalg.norm(ca_res1.get_coord() - ca_res2.get_coord())

def process_distances(df, structure):
    triple_links = []
    double_links = []

    triple_links_data = []
    double_links_data = []

    for _, row in df.iterrows():
        xl_values = ['XL A', 'XL B', 'XL C']
        chain_values = ['chain ID A', 'chain ID B', 'chain ID C']

        if pd.notna(row['XL A']) and pd.notna(row['XL B']) and pd.notna(row['XL C']):
            # Triple-links case
            for i in range(3):
                for j in range(i + 1, 3):
                    xl_i = row[xl_values[i]]
                    chain_i = row[chain_values[i]]
                    xl_j = row[xl_values[j]]
                    chain_j = row[chain_values[j]]

                    if pd.notna(xl_i) and pd.notna(xl_j) and pd.notna(chain_i) and pd.notna(chain_j):
                        res_i = int(re.sub(r'\D', '', xl_i))
                        res_j = int(re.sub(r'\D', '', xl_j))
                        distance = calculate_ca_distance(structure, res_i, chain_i, res_j, chain_j)
                        
                        if not np.isnan(distance):
                            triple_links.append(distance)
                            triple_links_data.append({
                                'chain ID 1': chain_i,
                                'res 1': res_i,
                                'chain ID 2': chain_j,
                                'res 2': res_j,
                                'distance': distance
                            })

        elif pd.notna(row['XL A']) and pd.notna(row['XL B']):
            # Double-links case
            xl_i = row['XL A']
            chain_i = row['chain ID A']
            xl_j = row['XL B']
            chain_j = row['chain ID B']

            if pd.notna(xl_i) and pd.notna(xl_j) and pd.notna(chain_i) and pd.notna(chain_j):
                res_i = int(re.sub(r'\D', '', xl_i))
                res_j = int(re.sub(r'\D', '', xl_j))
                distance = calculate_ca_distance(structure, res_i, chain_i, res_j, chain_j)

                if not np.isnan(distance):
                    double_links.append(distance)
                    double_links_data.append({
                        'chain ID 1': chain_i,
                        'res 1': res_i,
                        'chain ID 2': chain_j,
                        'res 2': res_j,
                        'distance': distance
                    })

    # Convert lists to DataFrames
    triple_links_df = pd.DataFrame(triple_links_data)
    double_links_df = pd.DataFrame(double_links_data)

    return triple_links, double_links, triple_links_df, double_links_df

def plot_histograms(triple_links, double_links):
    # Define figure dimensions in inches based on the golden ratio
    fig_width_inch = 7.2  # Width in inches (183 mm)
    fig_height_inch = fig_width_inch / 1.618  # Height in inches based on the golden ratio

    plt.figure(figsize=(fig_width_inch, fig_height_inch))

    # Histogram for triple-links
    plt.subplot(1, 2, 1)
    sns.histplot(triple_links, bins=30, kde=False, color='blue')
    plt.axvline(x=35.0, color='red', linestyle='--', label='x=35.0')
    plt.xlabel(r'Distance (${\AA}$)', fontsize=14)
    plt.ylabel('Frequency', fontsize=14)
    plt.title('Triple-links: All A, B, C', fontsize=16)
    plt.legend(loc='upper right', fontsize=12)

    # Calculate total number of triple-links
    total_triple_links = len(triple_links)
    # Add text annotation for triple-links
    plt.text(0.98, 0.85, f'Total Triple-links: {total_triple_links/3}', 
            ha='right', va='top', fontsize=12, color='blue', 
            transform=plt.gca().transAxes)

    # Histogram for double-links
    plt.subplot(1, 2, 2)
    sns.histplot(double_links, bins=30, kde=False, color='green')
    plt.axvline(x=35.0, color='red', linestyle='--', label='x=35.0')
    plt.xlabel(r'Distance (${\AA}$)', fontsize=14)
    plt.ylabel('Frequency', fontsize=14)
    plt.title('Double-links: A and B', fontsize=16)
    plt.legend(loc='upper right', fontsize=12)

    # Calculate total number of double-links
    total_double_links = len(double_links)
    # Add text annotation for double-links
    plt.text(0.98, 0.85, f'Total Double-links: {total_double_links}', 
            ha='right', va='top', fontsize=12, color='green', 
            transform=plt.gca().transAxes)

    plt.tight_layout()
    plt.savefig('xls_satisfaction_base.pdf', format='pdf', dpi=300)
    plt.show()
    
def replace_chain_ids(df, mapping):
    """Replaces chain IDs in the DataFrame with their corresponding gene names."""
    chain_id_columns = ['chain ID 1', 'chain ID 2']
    
    for col in chain_id_columns:
        if col in df.columns:
            df[col] = df[col].map(lambda x: mapping.get(x, x))
        else:
            print(f"Warning: Column '{col}' not found in the DataFrame.")
    # Rename the columns chain ID 1 and chain ID 2 to protein1 and protein2
    df.rename(columns={'chain ID 1': 'protein1', 'chain ID 2': 'protein2'}, inplace=True)

    return df

def main():
    # Load and filter the CSV file into a DataFrame
    df = pd.read_csv('base-of-proteasome.csv')

    alphabets = ['v', 'w', 'x', 'y', 'z', '0', '1']
    filtered_df = filter_chain_ids(df, alphabets)
    filtered_df.to_csv('filtered_base-of-proteasome.csv', index=False)

    # Load the PDB file
    pdb_file = 'base_proteasome.pdb'
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('structure', pdb_file)

    # Process the distances
    triple_links, double_links, triple_links_df, double_links_df = process_distances(filtered_df, structure)
    
    # Calculate statistics
    num_less_than_35_triple = sum(d < 35 for d in triple_links)
    num_less_than_35_double = sum(d < 35 for d in double_links)
    total_triple = len(triple_links)
    total_double = len(double_links)

    print(f"Triple-links: Number of distances less than 35: {num_less_than_35_triple}")
    print(f"Triple-links: Total number of calculated distances: {total_triple}")
    print(f"Double-links: Number of distances less than 35: {num_less_than_35_double}")
    print(f"Double-links: Total number of calculated distances: {total_double}")

    # Plot histograms
    plot_histograms(triple_links, double_links)

    # Save the results to CSV files
    triple_links_df.to_csv('triple_links_distances.csv', index=False)
    double_links_df.to_csv('double_links_distances.csv', index=False)
    
    # Load the DataFrame from CSV
    df = pd.read_csv('triple_links_distances.csv')
    # Define the replacement mapping
    replacement_mapping = {
        'x': 'Rpt6',
        'y': 'Rpt3',
        'z': 'Rpt4',
        '1': 'Rpn2',
        '0': 'Rpt5',
        'w': 'Rpt2',
        'v': 'Rpt1'
    }
    # Replace the chain IDs with the gene names
    df = replace_chain_ids(df, replacement_mapping)
    # Save the modified DataFrame to a new CSV file
    # Lose the last column before saving
    df = df.iloc[:, :-1]
    df.to_csv('replaced_triple_links.csv', index=False)
    # Replace the chain IDs with the gene names for double links
    double_links_df = replace_chain_ids(double_links_df, replacement_mapping)
    # Save the modified DataFrame to a new CSV file
    # Lose the last column before saving
    double_links_df = double_links_df.iloc[:, :-1]
    double_links_df.to_csv('replaced_double_links.csv', index=False)
    # find similar rows
    df1 = pd.read_csv('replaced_triple_links.csv')
    df2 = pd.read_csv('replaced_double_links.csv')
    common_rows = pd.merge(df1,df2, on=['protein1', 'res 1', 'protein2', 'res 2'], how='inner')
    print("common rows between the datasets: ", common_rows)
    
    # Initialize results list
    results = []

    # Process DataFrame in chunks of 3 rows
    for start in range(0, len(df), 3):
        chunk = df.iloc[start:start+3]
    
        # Collect unique (protein, res) pairs and their occurrences
        for _, row in chunk.iterrows():
            results.append([row['protein1'], row['res 1'], f'p{(start // 3) + 1}', 1])
            results.append([row['protein2'], row['res 2'], f'p{(start // 3) + 1}', 1])

    # Create DataFrame from results
    result_df = pd.DataFrame(results, columns=['Protein1', 'Residue1', 'Protein2', 'Residue2'])

    # Remove duplicates
    result_df = result_df.drop_duplicates()

    # Print the DataFrame
    #print(result_df)
    result_df.to_csv('final_triple_xls.csv', index=False)
    
if __name__ == '__main__':
    main()
