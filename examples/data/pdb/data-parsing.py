#!/usr/bin/env python3.11

import pandas as pd
import re
#from Bio.PDB import *
from Bio import PDB
from Bio.PDB import *
from Bio import SeqIO
import MDAnalysis as mda
import sys
import warnings

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
#processed_data['Score A'] = psome_df['Score A']
#processed_data['Score B'] = psome_df['Score B']
#processed_data['Score C'] = psome_df['Score C']
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
print(processed_data)
#-----------------------------------------------------------------------------
# Function to calculate the distance between two C-alpha atoms in a CIF file
#-----------------------------------------------------------------------------
def calculate_ca_distance(cif_file, res1_num, chain1_id, res2_num, chain2_id):
    # Create a CIF parser
    parser = PDB.MMCIFParser()
    
    # Load the structure from the CIF file
    structure = parser.get_structure('structure', cif_file)
    
    # Select the C-alpha atoms for the specified residues and chains
    ca_res1 = None
    ca_res2 = None
    
    for model in structure:
        for chain in model:
            if chain.id == chain1_id:
                for residue in chain:
                    if residue.get_id()[1] == res1_num and 'CA' in residue:
                        ca_res1 = residue['CA']
            if chain.id == chain2_id:
                for residue in chain:
                    if residue.get_id()[1] == res2_num and 'CA' in residue:
                        ca_res2 = residue['CA']
    
    if ca_res1 is None or ca_res2 is None:
        print("Error: One or both specified C-alpha atoms were not found.")
        return
    
    # Calculate the distance between the C-alpha atoms
    coord1 = ca_res1.get_coord()
    coord2 = ca_res2.get_coord()
    distance = np.linalg.norm(coord1 - coord2)
    
    # Print the distance
    print(f"Distance between C-alpha atoms: {distance:.2f} Å")

cif_file = './5gjr.cif'
res1_num = 666
chain1_id = 'N'
res2_num = 807
chain2_id = 'N'



calculate_ca_distance(cif_file, res1_num, chain1_id, res2_num, chain2_id)


