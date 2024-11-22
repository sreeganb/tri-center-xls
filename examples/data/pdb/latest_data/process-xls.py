#!/usr/bin/env python3.11

import pandas as pd
import os

def clean_residue(residue):
    """
    Remove all 'K' prefixes from a residue string.
    """
    if pd.isna(residue):
        return residue
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
    df = pd.read_csv(file_path)

    # Display columns
    print("Columns in the DataFrame:", df.columns.tolist())

    # Rename columns
    df.rename(columns={
        'Protein A': 'Protein1',
        'XL A': 'Residue1',
        'Gene A': 'Gene1',
        'Protein B': 'Protein2',
        'XL B': 'Residue2',
        'Gene B': 'Gene2',
        'Protein C': 'Protein3',
        'XL C': 'Residue3',
        'Gene C': 'Gene3',
    }, inplace=True)

    # Replace NaN values with empty strings
    df = df.fillna('')

    # Identify triple and double links
    df['is_triple'] = df['Residue3'].astype(bool)

    # Split into triple_links and double_links
    triple_links = df[df['is_triple']].copy()
    double_links = df[~df['is_triple']].copy()

    # Select relevant columns including Gene columns
    double_links = double_links[['Gene1', 'Residue1',
                                 'Gene2', 'Residue2']]
    triple_links = triple_links[[ 'Gene1', 'Residue1',
                                  'Gene2', 'Residue2',
                                  'Gene3', 'Residue3']]

    # Apply the clean_residue function to residue columns
    residue_columns = ['Residue1', 'Residue2', 'Residue3']
    for col in residue_columns:
        if col in double_links.columns:
            double_links[col] = double_links[col].apply(clean_residue)
        if col in triple_links.columns:
            triple_links[col] = triple_links[col].apply(clean_residue)

    # Replace special characters in protein columns
    protein_columns = ['Gene1', 'Gene2', 'Gene3']
    for col in protein_columns:
        if col in double_links.columns:
            double_links[col] = double_links[col].apply(replace_special_characters)
        if col in triple_links.columns:
            triple_links[col] = triple_links[col].apply(replace_special_characters)

    return triple_links, double_links

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

    # Save the ambiguous and cleaned DataFrames to CSV files
    output_dir = 'proteasome_data'
    os.makedirs(output_dir, exist_ok=True)

    ambiguous_double.to_csv(os.path.join(output_dir, 'ambiguous_double_links.csv'), index=False)
    cleaned_double_links.to_csv(os.path.join(output_dir, 'cleaned_double_links.csv'), index=False)
    ambiguous_triple.to_csv(os.path.join(output_dir, 'ambiguous_triple_links.csv'), index=False)
    cleaned_triple_links.to_csv(os.path.join(output_dir, 'cleaned_triple_links.csv'), index=False)

    return ambiguous_double, ambiguous_triple, cleaned_double_links, cleaned_triple_links

def filter_and_save_gene_entries(links, link_type='double', is_ambiguous=False):
    # Define the set of genes to filter
    gene_set = {'Rpt1', 'Rpt2', 'Rpt3', 'Rpt4', 'Rpt5', 'Rpt6', 'Rpn2'}

    # Filter rows where all gene columns are in the gene_set
    gene_columns = [col for col in links.columns if 'Gene' in col]
    filtered_links = links[
        links[gene_columns].apply(lambda row: all(gene in gene_set for gene in row), axis=1)
    ]

    # Determine the output directory and filename
    output_dir = 'proteasome_data'
    os.makedirs(output_dir, exist_ok=True)

    if is_ambiguous:
        filename = f'ambiguous_{link_type}_links_filtered.csv'
    else:
        filename = f'cleaned_{link_type}_links_filtered.csv'

    # Save the filtered DataFrame to a CSV file
    filtered_links.to_csv(os.path.join(output_dir, filename), index=False)

    return filtered_links

def main():
    file_path = 'input/data.csv'  # Update with the path to your CSV file

    # Process the CSV file
    triple_links, double_links = process_csv(file_path)

    # Filter ambiguous links and save the cleaned data
    ambiguous_double, ambiguous_triple, cleaned_double_links, cleaned_triple_links = filter_ambiguous_links(
        double_links, triple_links
    )

    # Filter and save gene entries for cleaned and ambiguous links
    cleaned_double_filtered = filter_and_save_gene_entries(
        cleaned_double_links, link_type='double', is_ambiguous=False
    )
    cleaned_triple_filtered = filter_and_save_gene_entries(
        cleaned_triple_links, link_type='triple', is_ambiguous=False
    )
    ambiguous_double_filtered = filter_and_save_gene_entries(
        ambiguous_double, link_type='double', is_ambiguous=True
    )
    ambiguous_triple_filtered = filter_and_save_gene_entries(
        ambiguous_triple, link_type='triple', is_ambiguous=True
    )

if __name__ == '__main__':
    main()