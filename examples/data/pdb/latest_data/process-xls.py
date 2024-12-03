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

def replace_special_characters(gene):
    """
    Replace special characters α and β with alpha and beta.
    """
    if pd.isna(gene):
        return gene
    return gene.replace('α', 'alpha').replace('β', 'beta')

def process_csv(file_path):
    # Read the CSV file into a DataFrame
    df = pd.read_csv(file_path)

    # Display columns
    print("Columns in the DataFrame:", df.columns.tolist())

    # Rename columns directly without swapping
    df.rename(columns={
        'Gene A': 'Gene1',
        'XL A': 'Residue1',
        'Gene B': 'Gene2',
        'XL B': 'Residue2',
        'Gene C': 'Gene3',
        'XL C': 'Residue3',
    }, inplace=True)

    # Replace NaN values with empty strings
    df = df.fillna('')

    # Identify triple and double links
    df['is_triple'] = df['Residue3'].astype(bool)

    # Split into triple_links and double_links
    triple_links = df[df['is_triple']].copy()
    double_links = df[~df['is_triple']].copy()

    # Select relevant columns (only genes and residues)
    double_links = double_links[['Gene1', 'Residue1', 'Gene2', 'Residue2']]
    triple_links = triple_links[['Gene1', 'Residue1', 'Gene2', 'Residue2', 'Gene3', 'Residue3']]

    # Apply the clean_residue function to residue columns
    residue_columns = ['Residue1', 'Residue2', 'Residue3']
    for col in residue_columns:
        if col in double_links.columns:
            double_links[col] = double_links[col].apply(clean_residue)
        if col in triple_links.columns:
            triple_links[col] = triple_links[col].apply(clean_residue)

    # Replace special characters in gene columns
    gene_columns = ['Gene1', 'Gene2', 'Gene3']
    for col in gene_columns:
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

    # Remove duplicate rows and report the number of duplicates removed
    for df_name, df in [('cleaned_double_links', cleaned_double_links), ('cleaned_triple_links', cleaned_triple_links),
                        ('ambiguous_double', ambiguous_double), ('ambiguous_triple', ambiguous_triple)]:
        initial_count = len(df)
        df.drop_duplicates(inplace=True)
        duplicates_removed = initial_count - len(df)
        if duplicates_removed > 0:
            print(f"Removed {duplicates_removed} duplicate rows from {df_name}.")

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

    # Get gene and residue columns present in the DataFrame
    gene_columns = [col for col in links.columns if 'Gene' in col]
    residue_columns = [col for col in links.columns if 'Residue' in col]

    # Ensure that gene_columns and residue_columns are not empty
    if not gene_columns or not residue_columns:
        print("No Gene or Residue columns found in the DataFrame.")
        return pd.DataFrame()  # Return an empty DataFrame

    # Filter rows where all gene columns are in the gene_set
    filtered_links = links[
        links[gene_columns].apply(lambda row: all(gene in gene_set for gene in row), axis=1)
    ]

    # Remove duplicate rows and report the number of duplicates removed
    initial_count = len(filtered_links)
    filtered_links.drop_duplicates(inplace=True)
    duplicates_removed = initial_count - len(filtered_links)
    if duplicates_removed > 0:
        print(f"Removed {duplicates_removed} duplicate rows from filtered {link_type} links (ambiguous={is_ambiguous}).")

    # Determine the output directory and filename
    output_dir = 'proteasome_data'
    os.makedirs(output_dir, exist_ok=True)

    if is_ambiguous:
        filename = f'ambiguous_{link_type}_links_filtered.csv'
    else:
        filename = f'cleaned_{link_type}_links_filtered.csv'

    # Save the filtered DataFrame to a CSV file
    available_columns = [col for col in gene_columns + residue_columns if col in filtered_links.columns]
    filtered_links.to_csv(os.path.join(output_dir, filename), index=False, columns=available_columns)

    return filtered_links

def pair_triples(triple_links):
    """
    Convert each triple link into three paired links.

    For each row in triple_links with columns:
    'Gene1', 'Residue1', 'Gene2', 'Residue2', 'Gene3', 'Residue3'

    Create three rows:
    1. 'Gene1', 'Residue1', 'Gene2', 'Residue2'
    2. 'Gene1', 'Residue1', 'Gene3', 'Residue3'
    3. 'Gene2', 'Residue2', 'Gene3', 'Residue3'
    """
    paired_links_list = []

    for idx, row in triple_links.iterrows():
        # First pair: Gene1 with Gene2
        pair1 = {
            'Gene1': row['Gene1'],
            'Residue1': row['Residue1'],
            'Gene2': row['Gene2'],
            'Residue2': row['Residue2']
        }
        # Second pair: Gene1 with Gene3
        pair2 = {
            'Gene1': row['Gene1'],
            'Residue1': row['Residue1'],
            'Gene2': row['Gene3'],
            'Residue2': row['Residue3']
        }
        # Third pair: Gene2 with Gene3
        pair3 = {
            'Gene1': row['Gene2'],
            'Residue1': row['Residue2'],
            'Gene2': row['Gene3'],
            'Residue2': row['Residue3']
        }
        paired_links_list.extend([pair1, pair2, pair3])

    # Create a DataFrame from the list of paired links
    paired_triple_links = pd.DataFrame(paired_links_list)

    # Remove duplicate rows and report the number of duplicates removed
    initial_count = len(paired_triple_links)
    paired_triple_links.drop_duplicates(inplace=True)
    duplicates_removed = initial_count - len(paired_triple_links)
    if duplicates_removed > 0:
        print(f"Removed {duplicates_removed} duplicate rows from paired triple links.")

    return paired_triple_links

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

    # Check if cleaned_triple_filtered is not empty before pairing
    if not cleaned_triple_filtered.empty:
        # Pair the cleaned triple links into double links
        paired_triple_links = pair_triples(cleaned_triple_filtered)

        # Save the paired triple links to a CSV file
        output_dir = 'proteasome_data'
        os.makedirs(output_dir, exist_ok=True)
        paired_triple_links.to_csv(os.path.join(output_dir, 'paired_triple_links.csv'), index=False)
    else:
        print("No cleaned triple links available for pairing.")

if __name__ == '__main__':
    main()