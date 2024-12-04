#!/usr/bin/env python3.11

import pandas as pd
import os

def clean_residue(residue):
    """
    Remove 'K' prefixes and strip spaces.
    """
    if pd.isna(residue) or residue == '':
        return ''
    # Remove 'K' prefix and strip spaces
    residue = residue.replace('K', '').strip()
    return residue

def replace_special_characters(protein):
    """
    Replace special characters α and β with 'alpha' and 'beta'.
    """
    if pd.isna(protein):
        return protein
    return protein.replace('α', 'alpha').replace('β', 'beta')

def explode_dataframe(df, columns, separators):
    """
    Explodes the DataFrame by splitting the specified columns on the given separators.

    Args:
        df (pd.DataFrame): The DataFrame to explode.
        columns (list): List of columns to split.
        separators (str): String of separators to split on.

    Returns:
        pd.DataFrame: The exploded DataFrame.
    """
    for column in columns:
        df = df.assign(**{column: df[column].str.split(separators)}).explode(column)
    return df

def process_csv(file_path):
    # Read the CSV file into a DataFrame
    df = pd.read_csv(file_path)

    print("Total rows in original dataframe:", len(df))

    # Display columns
    print("Columns in the DataFrame:", df.columns.tolist())

    # Rename columns directly
    df.rename(columns={
        'Protein A': 'Protein1',
        'XL A': 'Residue1',
        'Protein B': 'Protein2',
        'XL B': 'Residue2',
        'Protein C': 'Protein3',
        'XL C': 'Residue3',
    }, inplace=True)

    # Replace NaN values with empty strings
    df = df.fillna('')

    # Remove 'K' from residues and strip spaces
    residue_columns = ['Residue1', 'Residue2', 'Residue3']
    for col in residue_columns:
        df[col] = df[col].apply(clean_residue)

    # Replace special characters in protein columns
    protein_columns = ['Protein1', 'Protein2', 'Protein3']
    for col in protein_columns:
        df[col] = df[col].apply(replace_special_characters)

    return df

def save_dataframes(df, output_folder):
    # Create output directory if it doesn't exist
    os.makedirs(output_folder, exist_ok=True)

    total_rows = len(df)
    print(f"Total rows in dataframe: {total_rows}")

    # Identify triple and double links without adding 'is_triple' column
    triple_links = df[df['Residue3'].astype(bool)].copy()
    double_links = df[~df['Residue3'].astype(bool)].copy()

    print(f"Total triple links: {len(triple_links)}")
    print(f"Total double links: {len(double_links)}")

    # Remove 'Residue3' and 'Protein3' from double links
    double_links = double_links.drop(columns=['Protein3', 'Residue3'])

    # Save ambiguous and unique links separately
    # Define separators
    separators = r';|\|'  # split on ';' or '|'

    # Identify ambiguous rows (those with separators in any residue)
    residue_columns = ['Residue1', 'Residue2']
    triple_residue_columns = ['Residue1', 'Residue2', 'Residue3']

    # For double links
    ambiguous_double = double_links[
        double_links[residue_columns].apply(lambda x: x.str.contains(separators).any(), axis=1)
    ]
    unique_double = double_links.drop(ambiguous_double.index)

    print(f"Unique double links: {len(unique_double)}")
    print(f"Ambiguous double links: {len(ambiguous_double)}")
    print(f"Dropped rows in double links: {len(double_links) - len(unique_double) - len(ambiguous_double)}")

    # For triple links
    ambiguous_triple = triple_links[
        triple_links[triple_residue_columns].apply(lambda x: x.str.contains(separators).any(), axis=1)
    ]
    unique_triple = triple_links.drop(ambiguous_triple.index)

    print(f"Unique triple links: {len(unique_triple)}")
    print(f"Ambiguous triple links: {len(ambiguous_triple)}")
    print(f"Dropped rows in triple links: {len(triple_links) - len(unique_triple) - len(ambiguous_triple)}")

    # Save ambiguous and unique links
    # Ensure columns are in the specified order and unnecessary columns are removed
    double_columns_order = ['Protein1', 'Residue1', 'Protein2', 'Residue2']
    triple_columns_order = ['Protein1', 'Residue1', 'Protein2', 'Residue2', 'Protein3', 'Residue3']

    ambiguous_double.to_csv(
        os.path.join(output_folder, 'ambiguous_double_links.csv'),
        index=False,
        columns=double_columns_order
    )
    unique_double.to_csv(
        os.path.join(output_folder, 'unique_double_links.csv'),
        index=False,
        columns=double_columns_order
    )
    ambiguous_triple.to_csv(
        os.path.join(output_folder, 'ambiguous_triple_links.csv'),
        index=False,
        columns=triple_columns_order
    )
    unique_triple.to_csv(
        os.path.join(output_folder, 'unique_triple_links.csv'),
        index=False,
        columns=triple_columns_order
    )

    # Explode ambiguous links
    # For double links
    ambiguous_double_exploded = explode_dataframe(ambiguous_double, ['Residue1', 'Residue2'], separators)
    unique_double_exploded = unique_double.copy()

    # Combine unique and exploded ambiguous double links
    combined_double_links = pd.concat([unique_double_exploded, ambiguous_double_exploded], ignore_index=True)

    print(f"Total combined double links before dropping duplicates: {len(combined_double_links)}")
    combined_double_links.drop_duplicates(inplace=True)
    print(f"Total combined double links after dropping duplicates: {len(combined_double_links)}")

    # For triple links
    ambiguous_triple_exploded = explode_dataframe(ambiguous_triple, ['Residue1', 'Residue2', 'Residue3'], separators)
    unique_triple_exploded = unique_triple.copy()

    # Combine unique and exploded ambiguous triple links
    combined_triple_links = pd.concat([unique_triple_exploded, ambiguous_triple_exploded], ignore_index=True)

    print(f"Total combined triple links before dropping duplicates: {len(combined_triple_links)}")
    combined_triple_links.drop_duplicates(inplace=True)
    print(f"Total combined triple links after dropping duplicates: {len(combined_triple_links)}")

    # Save combined links
    combined_double_links.to_csv(
        os.path.join(output_folder, 'combined_double_links.csv'),
        index=False,
        columns=double_columns_order
    )
    combined_triple_links.to_csv(
        os.path.join(output_folder, 'combined_triple_links.csv'),
        index=False,
        columns=triple_columns_order
    )

    # Pair triples into doubles
    paired_triple_links = pair_triples(combined_triple_links)
    print(f"Total paired triple links (converted to doubles): {len(paired_triple_links)}")

    # Save paired triple links separately (do not combine with double links)
    paired_triple_links.to_csv(
        os.path.join(output_folder, 'paired_triple_links.csv'),
        index=False,
        columns=double_columns_order
    )

def pair_triples(triple_links):
    """
    Convert each triple link into three paired links.

    For each row in triple_links with columns:
    'Protein1', 'Residue1', 'Protein2', 'Residue2', 'Protein3', 'Residue3'

    Create three rows:
    1. 'Protein1', 'Residue1', 'Protein2', 'Residue2'
    2. 'Protein1', 'Residue1', 'Protein3', 'Residue3'
    3. 'Protein2', 'Residue2', 'Protein3', 'Residue3'
    """
    paired_links_list = []

    for idx, row in triple_links.iterrows():
        # First pair: Protein1 with Protein2
        pair1 = {
            'Protein1': row['Protein1'],
            'Residue1': row['Residue1'],
            'Protein2': row['Protein2'],
            'Residue2': row['Residue2']
        }
        # Second pair: Protein1 with Protein3
        pair2 = {
            'Protein1': row['Protein1'],
            'Residue1': row['Residue1'],
            'Protein2': row['Protein3'],
            'Residue2': row['Residue3']
        }
        # Third pair: Protein2 with Protein3
        pair3 = {
            'Protein1': row['Protein2'],
            'Residue1': row['Residue2'],
            'Protein2': row['Protein3'],
            'Residue2': row['Residue3']
        }
        paired_links_list.extend([pair1, pair2, pair3])

    # Create a DataFrame from the list of paired links
    paired_triple_links = pd.DataFrame(paired_links_list)

    # Remove duplicate rows and report the number of duplicates removed
    initial_count = len(paired_triple_links)
    paired_triple_links.drop_duplicates(inplace=True)
    duplicates_removed = initial_count - len(paired_triple_links)
    print(f"Removed {duplicates_removed} duplicate rows from paired triple links.")

    return paired_triple_links

if __name__ == '__main__':
    file_path = 'input_data/crosslinks.csv'  # Replace with your actual CSV file path
    df = process_csv(file_path)

    # Save full data
    print("\nProcessing full data...")
    save_dataframes(df, 'output_data/full_data')

    # Define the set of proteins to filter
    protein_set = {'Rpt1', 'Rpt2', 'Rpt3', 'Rpt4', 'Rpt5', 'Rpt6', 'Rpn2'}

    # Filter for base of the proteasome
    protein_columns = ['Protein1', 'Protein2', 'Protein3']
    filtered_df = df[
        df[protein_columns].apply(
            lambda row: all((protein in protein_set or protein == '') for protein in row), axis=1)
    ].copy()

    print(f"\nTotal rows after filtering for base of the proteasome: {len(filtered_df)}")
    print(f"Total rows dropped: {len(df) - len(filtered_df)}")

    # Save filtered data
    print("\nProcessing base of the proteasome data...")
    save_dataframes(filtered_df, 'output_data/base_proteasome')