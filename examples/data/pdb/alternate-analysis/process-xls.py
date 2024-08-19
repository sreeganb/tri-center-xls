#!/usr/bin/env python3.11

import pandas as pd
import os

def clean_residue(residue):
    """
    Remove all 'K' prefixes from a residue string.
    """
    if pd.isna(residue):
        return residue
    # Replace all occurrences of 'K' with an empty string
    return residue.replace('K', '')

def replace_special_characters(protein):
    """
    Replace special characters α and β with alpha and beta.
    """
    if pd.isna(protein):
        return protein
    return protein.replace('α', 'alpha').replace('β', 'beta')

def process_csv(file_path):
    # Read the CSV file into a DataFrame
    df = pd.read_csv(file_path, header=None)
    
    # Apply the clean_residue function to each relevant column
    for column in [1, 3, 5]:  # Columns Residue1 (1), Residue2 (3), Residue3 (5)
        df[column] = df[column].apply(clean_residue)
    
    # Replace special characters in protein columns
    for column in [0, 2, 4]:  # Columns protein1 (0), protein2 (2), protein3 (4)
        df[column] = df[column].apply(replace_special_characters)
    
    # Replace NaN values with empty strings
    df = df.fillna('')
    
    # Create triple_links DataFrame with all six columns non-empty
    triple_links = df[(df[0] != '') & (df[1] != '') & (df[2] != '') & (df[3] != '') & (df[4] != '') & (df[5] != '')]
    
    # Create double_links DataFrame with the first four columns non-empty and the last two columns empty
    double_links = df[(df[0] != '') & (df[1] != '') & (df[2] != '') & (df[3] != '') & (df[4] == '') & (df[5] == '')]
    
    # Drop the last two columns from double_links
    double_links = double_links.drop(columns=[4, 5])
    
    # Add column names to double_links
    double_links.columns = ['Protein1', 'Residue1', 'Protein2', 'Residue2']
    
    # Create the output directory if it doesn't exist
    output_dir = 'proteasome_data'
    os.makedirs(output_dir, exist_ok=True)
    
    # Save the double_links DataFrame to a new CSV file
    double_links.to_csv(os.path.join(output_dir, 'double_links.csv'), index=False)

    return triple_links, double_links

def filter_ambiguous_links(double_links, triple_links):
    # Columns to check for ambiguous characters
    double_columns_to_check = ['Residue1', 'Residue2']
    triple_columns_to_check = [1, 3, 5]
    
    # Filter ambiguous rows
    ambiguous_double = double_links[double_links[double_columns_to_check].apply(lambda x: x.str.contains(r';|\|').any(), axis=1)]
    ambiguous_triple = triple_links[triple_links[triple_columns_to_check].apply(lambda x: x.str.contains(r';|\|').any(), axis=1)]
    
    # Remove ambiguous rows from the original DataFrames
    cleaned_double_links = double_links.drop(ambiguous_double.index)
    cleaned_triple_links = triple_links.drop(ambiguous_triple.index)
    
    # Ensure the column names are correct before saving
    ambiguous_double.columns = ['Protein1', 'Residue1', 'Protein2', 'Residue2']
    cleaned_double_links.columns = ['Protein1', 'Residue1', 'Protein2', 'Residue2']
    ambiguous_triple.columns = ['Protein1', 'Residue1', 'Protein2', 'Residue2', 'Protein3', 'Residue3']
    cleaned_triple_links.columns = ['Protein1', 'Residue1', 'Protein2', 'Residue2', 'Protein3', 'Residue3']
    
    # Save the ambiguous and cleaned DataFrames to CSV files
    output_dir = 'proteasome_data'
    ambiguous_double.to_csv(os.path.join(output_dir, 'ambiguous_double_links.csv'), index=False)
    cleaned_double_links.to_csv(os.path.join(output_dir, 'cleaned_double_links.csv'), index=False)
    ambiguous_triple.to_csv(os.path.join(output_dir, 'ambiguous_triple_links.csv'), index=False)
    cleaned_triple_links.to_csv(os.path.join(output_dir, 'cleaned_triple_links.csv'), index=False)
    
    return ambiguous_double, ambiguous_triple, cleaned_double_links, cleaned_triple_links

def filter_and_save_protein_entries(links, link_type='double', is_ambiguous=False):
    # Define the set of proteins to filter
    protein_set = {'Rpt1', 'Rpt2', 'Rpt3', 'Rpt4', 'Rpt5', 'Rpt6', 'Rpn2'}
    
    # Filter the DataFrame based on the link type
    if link_type == 'double':
        filtered_df = links[
            links['Protein1'].isin(protein_set) & 
            links['Protein2'].isin(protein_set)
        ]
    elif link_type == 'triple':
        filtered_df = links[
            links['Protein1'].isin(protein_set) & 
            links['Protein2'].isin(protein_set) & 
            links['Protein3'].isin(protein_set)
        ]
    
    # Determine the output file name based on the link type and whether the links are ambiguous
    output_dir = 'base_subcomplex_data'
    os.makedirs(output_dir, exist_ok=True)
    if link_type == 'double':
        file_name = 'filtered_ambiguous_double_links.csv' if is_ambiguous else 'filtered_cleaned_double_links.csv'
    elif link_type == 'triple':
        file_name = 'filtered_ambiguous_triple_links.csv' if is_ambiguous else 'filtered_cleaned_triple_links.csv'
    
    # Save the filtered DataFrame to a file
    filtered_df.to_csv(os.path.join(output_dir, file_name), index=False)

def save_unique_links():
    output_dir = 'base_subcomplex_data'
    
    # Read the filtered cleaned double links file
    double_links_path = os.path.join(output_dir, 'filtered_cleaned_double_links.csv')
    if os.path.exists(double_links_path):
        double_links_df = pd.read_csv(double_links_path)
        unique_double_links_df = double_links_df.drop_duplicates()
        unique_double_links_df.to_csv(os.path.join(output_dir, 'unique_double_links.csv'), index=False)
    
    # Read the filtered cleaned triple links file
    triple_links_path = os.path.join(output_dir, 'filtered_cleaned_triple_links.csv')
    if os.path.exists(triple_links_path):
        triple_links_df = pd.read_csv(triple_links_path)
        unique_triple_links_df = triple_links_df.drop_duplicates()
        unique_triple_links_df.to_csv(os.path.join(output_dir, 'unique_triple_links.csv'), index=False)

def main():
    # Specify the path to your CSV file
    file_path = 'proteasome-data.csv'
    triple_links, double_links = process_csv(file_path)

    # Filter ambiguous links
    ambiguous_double, ambiguous_triple, cleaned_double_links, cleaned_triple_links = filter_ambiguous_links(double_links, triple_links)

    # Print the resulting DataFrames
    print("\nAmbiguous Double Links DataFrame:")
    print(ambiguous_double)
    print("\nAmbiguous Triple Links DataFrame:")
    print(ambiguous_triple)
    print("\nCleaned Double Links DataFrame:")
    print(cleaned_double_links)
    print("\nCleaned Triple Links DataFrame:")
    print(cleaned_triple_links)
    
    # Filter and save the cleaned double links
    filter_and_save_protein_entries(cleaned_double_links, link_type='double', is_ambiguous=False)
    
    # Filter and save the ambiguous double links
    filter_and_save_protein_entries(ambiguous_double, link_type='double', is_ambiguous=True)
    
    # Filter and save the cleaned triple links
    filter_and_save_protein_entries(cleaned_triple_links, link_type='triple', is_ambiguous=False)
    
    # Filter and save the ambiguous triple links
    filter_and_save_protein_entries(ambiguous_triple, link_type='triple', is_ambiguous=True)
    
    # Save unique links
    save_unique_links()

if __name__ == "__main__":
    main()
