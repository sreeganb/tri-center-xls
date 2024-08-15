#!/usr/bin/env python3.12

import pandas as pd
import re
import os

def clean_residue(entry):
    """Removes the leading 'K' from numeric parts of a residue.

    Args:
        entry: The input string.

    Returns:
        The modified string.
    """
    # if entry contains | then print it out 
    parts = entry.split(';')
    for i in range(len(parts)):
        subparts = parts[i].split('|')
        parts[i] = '|'.join(subpart[1:] if subpart.startswith('K') else subpart for subpart in subparts)
    if ';' in entry:
        print(entry)
        print(';'.join(parts))
    
    return ';'.join(parts)

def parse_protein_residue(protein_residue_str):
    """
    Parses a protein residue string into separate protein and residue components.
    
    Args:
        protein_residue_str (str): The protein residue string to be parsed.
        
    Returns:
        list: A list containing three proteins and their corresponding residues.
    """
    if '-' in protein_residue_str:
        parts = re.split('[-|]', protein_residue_str)
        proteins_residues = [part.split(':') for part in parts]
    else:
        proteins_residues = [part.split(':') for part in protein_residue_str.split('|')]

    proteins = []
    residues = []

    for pr in proteins_residues:
        if len(pr) == 2:
            protein, residue = pr
            proteins.append(protein)
            residues.append(clean_residue(residue))
        else:
            proteins.append('')
            residues.append('')

    while len(proteins) < 3:
        proteins.append('')
    while len(residues) < 3:
        residues.append('')

    return [
        proteins[0], residues[0] if len(residues) > 0 else '',
        proteins[1], residues[1] if len(residues) > 1 else '',
        proteins[2], residues[2] if len(residues) > 2 else ''
    ]

def extract_protein_residue(row):
    """
    Extracts protein and residue information from a CSV row.
    
    Args:
        row (str): The CSV row to be processed.
        
    Returns:
        list: A list containing three proteins and their corresponding residues.
    """
    components = row.split(',')
    if len(components) < 5:
        return ['', '', '', '', '', '']

    proteins = [comp.split('|')[2] if len(comp.split('|')) >= 3 else '' for comp in components[1:4]]
    residues = [clean_residue(res) for res in components[4:] if res]
    
    while len(proteins) < 3:
        proteins.append('')
    while len(residues) < 3:
        residues.append('')

    return [
        proteins[0], residues[0],
        proteins[1], residues[1],
        proteins[2], residues[2]
    ]

def process_csv(file_path):
    """
    Processes a CSV file and replaces special characters in the data.
    
    Args:
        file_path (str): The path to the CSV file to be processed.
        
    Returns:
        tuple: A tuple containing the processed DataFrame, triple links DataFrame, double links DataFrame, ambiguous triples DataFrame, and ambiguous doubles DataFrame.
    """
    df_raw = pd.read_csv(file_path, header=None)
    
    # Replace special characters
    df_raw = df_raw.applymap(lambda x: x.replace('β', 'beta').replace('α', 'alpha') if isinstance(x, str) else x)
    
    # Apply extraction function to each row
    extracted_rows = []
    for index, row in df_raw.iterrows():
        row_data = row[0]
        print(row_data)
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

    # Filter out rows with ';' or '|' in Residue columns and add to ambiguous DataFrames
    ambiguous_triples = df_triple_links[
        df_triple_links[['Residue1', 'Residue2', 'Residue3']].apply(lambda x: x.str.contains(r';|\|').any(), axis=1)
    ]
    ambiguous_doubles = df_double_links[
        df_double_links[['Residue1', 'Residue2']].apply(lambda x: x.str.contains(r';|\|').any(), axis=1)
    ]

    df_triple_links = df_triple_links[~df_triple_links.index.isin(ambiguous_triples.index)]
    df_double_links = df_double_links[~df_double_links.index.isin(ambiguous_doubles.index)]

    # Ensure the directory exists
    output_dir = 'proteasome_data'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Save DataFrames to CSV files
    df_triple_links.to_csv(os.path.join(output_dir, 'triple-links.csv'), index=False)
    df_double_links.to_csv(os.path.join(output_dir, 'double-links.csv'), index=False)
    ambiguous_triples.to_csv(os.path.join(output_dir, 'ambiguous-triples.csv'), index=False)
    ambiguous_doubles.to_csv(os.path.join(output_dir, 'ambiguous-doubles.csv'), index=False)

    return df_processed, df_triple_links, df_double_links, ambiguous_triples, ambiguous_doubles

def filter_proteins(df_triple_links, df_double_links):
    """
    Filters the DataFrames to include only rows where Protein1, Protein2, and Protein3 (for triple links)
    or Protein1 and Protein2 (for double links) are exclusively from the allowed set.
    
    Args:
        df_triple_links (DataFrame): The DataFrame containing triple links.
        df_double_links (DataFrame): The DataFrame containing double links.
        
    Returns:
        tuple: A tuple containing the filtered triple links DataFrame and double links DataFrame.
    """
    allowed_proteins = {'Rpt1', 'Rpt2', 'Rpt3', 'Rpt4', 'Rpt5', 'Rpt6', 'Rpn2'}
    
    filtered_triple_links = df_triple_links[
        df_triple_links[['Protein1', 'Protein2', 'Protein3']].apply(lambda x: set(x).issubset(allowed_proteins), axis=1)
    ]
    
    filtered_double_links = df_double_links[
        df_double_links[['Protein1', 'Protein2']].apply(lambda x: set(x).issubset(allowed_proteins), axis=1)
    ]
    
    # Ensure the directory exists
    output_dir = 'base_subcomplex_data'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Save filtered DataFrames to CSV files
    filtered_triple_links.to_csv(os.path.join(output_dir, 'filtered-triple-links.csv'), index=False)
    filtered_double_links.to_csv(os.path.join(output_dir, 'filtered-double-links.csv'), index=False)
    
    return filtered_triple_links, filtered_double_links

def add_frequency_column(df):
    """
    Adds a frequency column to the DataFrame indicating the number of times each unique row appears.
    
    Args:
        df (DataFrame): The input DataFrame.
        
    Returns:
        DataFrame: A new DataFrame with unique rows and a frequency column.
    """
    # Group by all columns and count the frequency of each unique row
    df_with_frequency = df.groupby(list(df.columns)).size().reset_index(name='frequency')
    return df_with_frequency

def process_dataframes(df_triple_links, df_double_links):
    """
    Processes the given DataFrames to find unique rows and add a frequency column.
    
    Args:
        df_triple_links (DataFrame): The DataFrame containing triple links.
        df_double_links (DataFrame): The DataFrame containing double links.
        
    Returns:
        tuple: A tuple containing the unique triple links DataFrame and unique double links DataFrame with frequency columns.
    """
    # Find unique rows and add frequency column for df_triple_links
    unique_triple_links = add_frequency_column(df_triple_links)

    # Find unique rows and add frequency column for df_double_links
    unique_double_links = add_frequency_column(df_double_links)

    return unique_triple_links, unique_double_links

def main():
    """
    Main function to execute the script.
    """
    # Specify the path to your CSV file
    file_path = 'proteasome-data.csv'
    
    # Process the CSV file
    df_processed, df_triple_links, df_double_links, ambiguous_triples, ambiguous_doubles = process_csv(file_path)
    
    # Filter the DataFrames
    filtered_triple_links, filtered_double_links = filter_proteins(df_triple_links, df_double_links)
    # Example usage
    unique_triple_links, unique_double_links = process_dataframes(filtered_triple_links, filtered_double_links)

    # Print the resulting DataFrames
    print("Unique Triple Links DataFrame with Frequency:")
    print(unique_triple_links)
    print("\nUnique Double Links DataFrame with Frequency:")
    print(unique_double_links)
    
if __name__ == "__main__":
    main()