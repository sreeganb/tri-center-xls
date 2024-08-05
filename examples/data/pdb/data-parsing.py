#!/usr/bin/env python

import pandas as pd
import re
#from Bio.PDB import *
from Bio import PDB
from Bio.PDB import *
from Bio import SeqIO
import MDAnalysis as mda
import sys
import warnings
import numpy as np

# Suppress all warnings
warnings.filterwarnings("ignore")

# Extract chain ID and protein name information from FASTA file downloaded from RCSB pdb
def extract_protein_info_to_dataframe(fasta_file):
  """Extracts protein name and chain ID from a FASTA file into a DataFrame.

  Args:
    fasta_file: Path to the FASTA file.

  Returns:
    A pandas DataFrame with columns 'protein_name' and 'chain_id'.
  """

  data = []
  for record in SeqIO.parse(fasta_file, "fasta"):
    header = record.description
    # Here, you'll need to parse the header to extract protein name and chain ID
    # This is a simplified example, assuming a specific header format
    try:
      protein_name, chain_id = header.split("|")[2], header.split("|")[1]
    except IndexError:
      protein_name, chain_id = "Unknown", "Unknown"
    auth_pattern = r"\[auth (\w+)\]"
    chain_id = re.findall(auth_pattern, header)
    chain_id = ','.join(chain_id)
    data.append({'protein_name': protein_name, 'chain_id': chain_id})
    
  df = pd.DataFrame(data)
  return df

# Example usage:
fasta_file = "../fasta/rcsb_pdb_5GJR.fasta"
protein_df = extract_protein_info_to_dataframe(fasta_file)
print(protein_df)

protein_df.to_csv('output.csv', index=False)

def match_and_add_chain_ids(df, psome_df):
  """
  Matches proteins in df with those in psome_df and adds chain ID columns.

  Args:
    df: The DataFrame containing protein and chain ID information.
    psome_df: The DataFrame containing protein names.
  """

  psome_df['chain ID A'] = ''
  psome_df['chain ID B'] = ''
  psome_df['chain ID C'] = ''

  for index, row in psome_df.iterrows():
    protein_a = str(row['Subunit Name A']).lower()
    protein_b = str(row['Subunit Name B']).lower()
    protein_c = str(row['Subunit Name C']).lower()

    for _, df_row in df.iterrows():
        if df_row['protein_name'].lower() == protein_a.lower():
            psome_df.at[index, 'chain ID A'] = df_row['chain_id']
        if df_row['protein_name'].lower() == protein_b.lower():
            psome_df.at[index, 'chain ID B'] = df_row['chain_id']
        if df_row['protein_name'].lower() == protein_c.lower():
            psome_df.at[index, 'chain ID C'] = df_row['chain_id']
  return psome_df

df = pd.read_csv('output.csv')
psome_file = '../../derived_data/xl/proteasome-tri-crosslinks.csv'
psome_df = pd.read_csv(psome_file)
psome_df = match_and_add_chain_ids(df, psome_df)
processed_data = psome_df.iloc[:, :12]  # Select the first 12 columns
processed_data['Average Score'] = psome_df[['Score A', 'Score B', 'Score C']].mean(axis=1)
processed_data['chain ID A'] = psome_df['chain ID A']
processed_data['chain ID B'] = psome_df['chain ID B']
processed_data['chain ID C'] = psome_df['chain ID C']

def remove_rows_with_empty_chain_id_a(df):
  """Removes rows where 'chain ID A' is empty.

  Args:
    df: The DataFrame to process.

  Returns:
    The DataFrame with rows removed.
  """

  df = df.copy()
  df = df[df['chain ID A'] != '']
  return df

# Assuming processed_data is your DataFrame
processed_data = remove_rows_with_empty_chain_id_a(processed_data)
processed_data.to_csv('processed-data.csv', index=False)
# Further processing of the data
import csv

def remove_duplicate_rows(input_file, output_file):
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

# Example usage
input_file = 'processed-data.csv'
output_file = 'unique-xls.csv'
remove_duplicate_rows(input_file, output_file)
df1 = pd.read_csv('unique-xls.csv')
# Split the 'XL A' and 'XL B' columns by ';' and expand into separate rows
df1 = df1.assign(**{
    'XL A': df1['XL A'].str.split(r'[;|]'),
    'XL B': df1['XL B'].str.split(r'[;|]'),
    'XL C': df1['XL C'].str.split(r'[;|]')
}).explode('XL A').explode('XL B').explode('XL C')
#df1['XL A'] = df1['XL A'].str.replace('K', '', regex=False)
#df1['XL B'] = df1['XL B'].str.replace('K', '', regex=False)
#df1['XL C'] = df1['XL C'].str.replace('K', '', regex=False)
df1.to_csv('unique-xls.csv', index=False)
#-----------------------------------------------------------------------------
# Function to calculate the distance between two C-alpha atoms in a CIF file
#-----------------------------------------------------------------------------
def calculate_ca_distance(structure, res1_num, chain1_id, res2_num, chain2_id):
    # Select the C-alpha atoms for the specified residues and chains
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

    # Calculate the distance between the C-alpha atoms
    coord1 = ca_res1.get_coord()
    coord2 = ca_res2.get_coord()
    distance = np.linalg.norm(coord1 - coord2)

    # Print the distance
    #print(f"Distance between C-alpha atoms: {distance:.2f} Å")
    return distance

# Read the DataFrame
df1 = pd.read_csv('unique-xls.csv')
cif_file = './5gjr.cif'
# Create a CIF parser
parser = PDB.MMCIFParser()

# Load the structure from the CIF file
structure = parser.get_structure('structure', cif_file)
calpha_distances = []
non_empty_entries = df1['XL C'].notnull().sum()
print(f"Number of non-empty entries: {non_empty_entries}")
not_found = 0
for i in range(non_empty_entries):
  print(df1.loc[i, 'XL A'])
  if not df1.loc[i, 'XL A'].startswith('K') or not df1.loc[i, 'XL B'].startswith('K') or not df1.loc[i, 'XL C'].startswith('K'):
    continue
  res1_num = int(df1.loc[i, 'XL A'][1:])  # Drop the first character and convert to integer
  res2_num = int(df1.loc[i, 'XL B'][1:])  # Drop the first character and convert to integer
  res3_num = int(df1.loc[i, 'XL C'][1:])  # Drop the first character and convert to integer
  chain1_id = df1.loc[i, 'chain ID A'].split(',')[0]
  chain2_id = df1.loc[i, 'chain ID B'].split(',')[0]
  chain3_id = df1.loc[i, 'chain ID C'].split(',')[0]
  dis = calculate_ca_distance(structure, res1_num, chain1_id, res2_num, chain2_id)
  if dis != -1:
    calpha_distances.append(dis)
  else:
    not_found += 1
  dis = calculate_ca_distance(structure, res1_num, chain1_id, res3_num, chain3_id)
  if dis != -1:
    calpha_distances.append(dis)
  else:
    not_found += 1
  dis = calculate_ca_distance(structure, res2_num, chain2_id, res3_num, chain3_id)
  if dis != -1:
    calpha_distances.append(dis)
  else:
    not_found += 1
print(f"Number of not found entries: {not_found}")

# Plot a histogram of the C-alpha distances
import matplotlib.pyplot as plt

# Plot the histogram using Matplotlib

plt.hist(calpha_distances, bins=50, edgecolor='black')
plt.axvline(x=35.0, color='red', linestyle='--')  # Add a vertical line at x = 35.0
plt.xlabel('C-alpha Distance (Å)')
plt.ylabel('Frequency')
plt.title(r"Distribution C$\alpha$ distances of trilinks")
plt.tight_layout()  # Add this line to improve spacing between subplots
plt.savefig('triple-links.pdf')
plt.show()