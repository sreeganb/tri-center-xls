#!/usr/bin/env python

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

def clean_residue(residue):
    """
    Remove all 'K' prefixes from a residue string.
    """
    if pd.isna(residue):
        return residue
    return residue.replace('K', '')

def filter_ambiguous_links(double_links, triple_links):
    # Columns to check for ambiguous characters
    double_columns_to_check = ['Residue1', 'Residue2']
    triple_columns_to_check = ['Residue1', 'Residue2', 'Residue3']

    # Filter ambiguous rows
    ambiguous_double = double_links[
        double_links[double_columns_to_check].apply(lambda x: x.str.contains(r';|\|').any(), axis=1)
    ]
    ambiguous_triple = triple_links[
        triple_links[triple_columns_to_check].apply(lambda x: x.str.contains(r';|\|').any(), axis=1)
    ]

    # Remove ambiguous rows from the original DataFrames
    cleaned_double_links = double_links.drop(ambiguous_double.index)
    cleaned_triple_links = triple_links.drop(ambiguous_triple.index)

    # Remove duplicate rows and report the number of duplicates removed
    for df_name, df in [('cleaned_double_links', cleaned_double_links), ('cleaned_triple_links', cleaned_triple_links),
                        ('ambiguous_double', ambiguous_double), ('ambiguous_triple', ambiguous_triple)]:
        initial_count = len(df)
        df.drop_duplicates(inplace=True)
        duplicates_removed = initial_count - len(df)
        if duplicates_removed > 0:
            print(f"Removed {duplicates_removed} duplicate rows from {df_name}.")

    # Save the ambiguous and cleaned DataFrames to CSV files
    output_dir = 'output_data'
    os.makedirs(output_dir, exist_ok=True)

    ambiguous_double.to_csv(os.path.join(output_dir, 'ambiguous_double_links.csv'), index=False)
    cleaned_double_links.to_csv(os.path.join(output_dir, 'cleaned_double_links.csv'), index=False)
    ambiguous_triple.to_csv(os.path.join(output_dir, 'ambiguous_triple_links.csv'), index=False)
    cleaned_triple_links.to_csv(os.path.join(output_dir, 'cleaned_triple_links.csv'), index=False)

    return ambiguous_double, ambiguous_triple, cleaned_double_links, cleaned_triple_links

# Read in the data, it is in the form Protein1,Residue1,Protein2,Residue2,Protein3,Residue3
df_raw = pd.read_csv('input_data/csn_xls_raw.csv')

# Filter out rows that contain the entry SERBP1 in either Protein1, Protein2, or Protein3
# Replace NaN with empty string
df_raw = df_raw.fillna('')
df_filtered = df_raw[(df_raw['Protein1'] == 'SERBP1') | (df_raw['Protein2'] == 'SERBP1') | (df_raw['Protein3'] == 'SERBP1')]

# Save to a CSV file, make sure that the output directory exists
os.makedirs('output_data', exist_ok=True)
# df_filtered.to_csv('output_data/serbp1_interactions.csv', index=False)

# Separate the data into two, one where Protein3 and Residue3 are empty and one where Protein3 and Residue3 are not empty
# They will be called bifunctional and trifunctional respectively
df_bifunctional = df_filtered[(df_filtered['Protein3'] == '') & (df_filtered['Residue3'] == '')].copy()
df_trifunctional = df_filtered[(df_filtered['Protein3'] != '') & (df_filtered['Residue3'] != '')].copy()

# Save the data and drop repeated rows
df_bifunctional.drop_duplicates(inplace=True)
df_trifunctional.drop_duplicates(inplace=True)

# Drop columns Protein3 and Residue3 from the bifunctional data
df_bifunctional.drop(columns=['Protein3', 'Residue3'], inplace=True)

# Apply the clean_residue function to residue columns
residue_columns = ['Residue1', 'Residue2', 'Residue3']
for col in residue_columns:
    if col in df_bifunctional.columns:
        df_bifunctional[col] = df_bifunctional[col].apply(clean_residue)
    if col in df_trifunctional.columns:
        df_trifunctional[col] = df_trifunctional[col].apply(clean_residue)

filter_ambiguous_links(df_bifunctional, df_trifunctional)

# Save the bifunctional and trifunctional data to CSV files
df_bifunctional.to_csv('output_data/serbp1_bifunctional.csv', index=False)
df_trifunctional.to_csv('output_data/serbp1_trifunctional.csv', index=False)