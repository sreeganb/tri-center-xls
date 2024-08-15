#!/usr/bin/env python3.12

import pandas as pd
import re

def clean_residue(residue):
    # Split by semicolon if there are multiple residues
    segments = residue.split(';')
    # Clean each segment by removing the 'K' prefix
    cleaned_segments = [seg.lstrip('K') for seg in segments]
    # Join the cleaned segments back with semicolons
    return ';'.join(cleaned_segments)

def parse_protein_residue(protein_residue_str):
    # Handle cases where the format is complex
    if '-' in protein_residue_str:
        parts = re.split('[-|]', protein_residue_str)
        proteins_residues = [part.split(':') for part in parts]
    else:
        proteins_residues = [part.split(':') for part in protein_residue_str.split('|')]

    # Extract proteins and residues
    proteins = []
    residues = []

    for pr in proteins_residues:
        if len(pr) == 2:
            protein, residue = pr
            proteins.append(protein)
            residues.append(clean_residue(residue))
        else:
            # Handle cases with missing residue information
            proteins.append('')
            residues.append('')

    # Ensure proteins and residues lists have 3 elements each
    while len(proteins) < 3:
        proteins.append('')
    while len(residues) < 3:
        residues.append('')
    
    return [
        proteins[0], residues[0],
        proteins[1], residues[1],
        proteins[2], residues[2]
    ]

def extract_protein_residue(row):
    components = row.split(',')
    # Handle the case where there are fewer components
    if len(components) < 5:
        return ['', '', '', '', '', '']

    # Parse proteins and residues
    proteins = [comp.split('|')[2] if len(comp.split('|')) >= 3 else '' for comp in components[1:4]]
    residues = [clean_residue(res) for res in components[4:] if res]  # Remove empty residues and clean them
    
    # Ensure proteins and residues lists have 3 elements each
    while len(proteins) < 3:
        proteins.append('')
    while len(residues) < 3:
        residues.append('')

    # Combine proteins and residues into the final output
    return [
        proteins[0], residues[0] if len(residues) > 0 else '',
        proteins[1], residues[1] if len(residues) > 1 else '',
        proteins[2], residues[2] if len(residues) > 2 else ''
    ]

def process_csv(file_path):
    df_raw = pd.read_csv(file_path, header=None)
    
    # Apply extraction function to each row
    extracted_rows = []
    for index, row in df_raw.iterrows():
        row_data = row[0]
        if '-' in row_data or '|' in row_data:
            extracted_rows.append(parse_protein_residue(row_data))
        else:
            extracted_rows.append(extract_protein_residue(row_data))
    
    # Create a DataFrame with the desired columns
    df_processed = pd.DataFrame(extracted_rows, columns=['Protein1', 'Residue1', 'Protein2', 'Residue2', 'Protein3', 'Residue3'])
    
    # Split into two DataFrames based on the presence of non-empty entries
    df_triple_links = df_processed.dropna().loc[
        (df_processed[['Protein1', 'Residue1', 'Protein2', 'Residue2', 'Protein3', 'Residue3']] != '').all(axis=1)
    ]
    df_double_links = df_processed[~df_processed.index.isin(df_triple_links.index)]
    df_double_links = df_double_links.drop(columns=['Protein3', 'Residue3'])

    # Save DataFrames to CSV files
    df_triple_links.to_csv('triple-links.csv', index=False)
    df_double_links.to_csv('double-links.csv', index=False)

    return df_processed

# Specify the path to your CSV file
file_path = 'proteasome-data.csv'
df_processed = process_csv(file_path)

# Print the resulting DataFrames
print("Triple Links DataFrame:")
print(df_processed[(df_processed[['Protein1', 'Residue1', 'Protein2', 'Residue2', 'Protein3', 'Residue3']] != '').all(axis=1)])
print("\nDouble Links DataFrame:")
print(df_processed[~df_processed.index.isin(df_processed[(df_processed[['Protein1', 'Residue1', 'Protein2', 'Residue2', 'Protein3', 'Residue3']] != '').all(axis=1)].index)].drop(columns=['Protein3', 'Residue3']))
