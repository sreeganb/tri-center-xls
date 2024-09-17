import pandas as pd
import os

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

def filter_proteins(df_triple_xls, df_double_xls):                                                                   
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
                                                                                                                         
    filtered_triple_links = df_triple_xls[                                                                             
        df_triple_xls[['Protein1', 'Protein2', 'Protein3']].apply(lambda x: set(x).issubset(allowed_proteins), axis=1) 
    ]                                                                                                                    
                                                                                                                         
    filtered_double_links = df_double_xls[                                                                             
        df_double_xls[['Protein1', 'Protein2']].apply(lambda x: set(x).issubset(allowed_proteins), axis=1)             
    ]                                                                                                                    
                                                                                                                         
    return filtered_triple_links, filtered_double_links                      

def analyze_repeating_rows(file_path):
    """
    Analyzes a CSV file to identify repeating rows based on specified criteria.

    Args:
        file_path (str): The path to the CSV file.

    Returns:
        None
    """

    # Read the CSV file into a pandas DataFrame
    df = pd.read_csv(file_path)

    # Replace special characters in 'Common/Gene Name A', 'Common/Gene Name B', and 'Common/Gene Name C'
    df['Common/Gene Name A'] = df['Common/Gene Name A'].str.replace('α', 'alpha').str.replace('β', 'beta')
    df['Common/Gene Name B'] = df['Common/Gene Name B'].str.replace('α', 'alpha').str.replace('β', 'beta')
    df['Common/Gene Name C'] = df['Common/Gene Name C'].str.replace('α', 'alpha').str.replace('β', 'beta')

    # Handle missing values in 'XL C'
    df['XL C'] = df['XL C'].fillna('')

    # Create a column to indicate whether 'XL A', 'XL B', or 'XL C' contains ';' or '|'
    df['Is_Ambiguous'] = df['XL A'].str.contains(r';|\|') | df['XL B'].str.contains(r';|\|') | df['XL C'].str.contains(r';|\|')

    # Group by desired columns and count occurrences
    grouped_data = df.groupby(['Common/Gene Name A', 'XL A', 'Common/Gene Name B', 'XL B', 'Common/Gene Name C', 'XL C']).size().reset_index(name='Count')

    # Filter rows where the count is greater than 1
    repeating_rows = grouped_data[grouped_data['Count'] > 0]

    # Count the number of unique rows and repeating rows
    num_repeating_rows = len(repeating_rows)
    num_unique_rows = df[['XL A', 'XL B', 'XL C']].drop_duplicates().shape[0]

    # Output the results
    print("\nRepeating Rows:")
    print(repeating_rows.to_markdown(index=False))

    print(f"\nTotal Repeating Rows: {repeating_rows['Count'].sum()}")
    print(f"Number of Unique Rows: {num_unique_rows}")

    # Separate data based on ambiguity and 'XL C'
    ambiguous_triple_df = df[df['Is_Ambiguous'] & (df['XL C'] != '')]
    ambiguous_double_df = df[df['Is_Ambiguous'] & (df['XL C'] == '')]
    non_ambiguous_triple_df = df[~df['Is_Ambiguous'] & (df['XL C'] != '')]
    non_ambiguous_double_df = df[~df['Is_Ambiguous'] & (df['XL C'] == '')]

    # Explode the ambiguous DataFrames
    ambiguous_triple_df = explode_dataframe(ambiguous_triple_df, ['XL A', 'XL B', 'XL C'], r';|\|')
    ambiguous_double_df = explode_dataframe(ambiguous_double_df, ['XL A', 'XL B'], r';|\|')

    # Group and count for each subset
    grouped_ambiguous_triple_data = ambiguous_triple_df.groupby(['Common/Gene Name A', 'XL A', 'Common/Gene Name B', 'XL B', 'Common/Gene Name C', 'XL C']).size().reset_index(name='Count')
    grouped_ambiguous_double_data = ambiguous_double_df.groupby(['Common/Gene Name A', 'XL A', 'Common/Gene Name B', 'XL B']).size().reset_index(name='Count')
    grouped_non_ambiguous_triple_data = non_ambiguous_triple_df.groupby(['Common/Gene Name A', 'XL A', 'Common/Gene Name B', 'XL B', 'Common/Gene Name C', 'XL C']).size().reset_index(name='Count')
    grouped_non_ambiguous_double_data = non_ambiguous_double_df.groupby(['Common/Gene Name A', 'XL A', 'Common/Gene Name B', 'XL B']).size().reset_index(name='Count')

    # Filter repeating rows for each subset
    repeating_ambiguous_triple_rows = grouped_ambiguous_triple_data[grouped_ambiguous_triple_data['Count'] > 0]
    repeating_ambiguous_double_rows = grouped_ambiguous_double_data[grouped_ambiguous_double_data['Count'] > 0]
    repeating_non_ambiguous_triple_rows = grouped_non_ambiguous_triple_data[grouped_non_ambiguous_triple_data['Count'] > 0]
    repeating_non_ambiguous_double_rows = grouped_non_ambiguous_double_data[grouped_non_ambiguous_double_data['Count'] > 0]

    # Count unique rows for each subset
    num_repeating_ambiguous_triple_rows = len(repeating_ambiguous_triple_rows)
    num_unique_ambiguous_triple_rows = ambiguous_triple_df[['XL A', 'XL B', 'XL C']].drop_duplicates().shape[0]

    num_repeating_ambiguous_double_rows = len(repeating_ambiguous_double_rows)
    num_unique_ambiguous_double_rows = ambiguous_double_df[['XL A', 'XL B', 'XL C']].drop_duplicates().shape[0]

    num_repeating_non_ambiguous_triple_rows = len(repeating_non_ambiguous_triple_rows)
    num_unique_non_ambiguous_triple_rows = non_ambiguous_triple_df[['XL A', 'XL B', 'XL C']].drop_duplicates().shape[0]

    num_repeating_non_ambiguous_double_rows = len(repeating_non_ambiguous_double_rows)
    num_unique_non_ambiguous_double_rows = non_ambiguous_double_df[['XL A', 'XL B', 'XL C']].drop_duplicates().shape[0]

    # Calculate total counts for each subset
    total_repeating_ambiguous_triple_rows = repeating_ambiguous_triple_rows['Count'].sum()
    total_repeating_ambiguous_double_rows = repeating_ambiguous_double_rows['Count'].sum()
    total_repeating_non_ambiguous_triple_rows = repeating_non_ambiguous_triple_rows['Count'].sum()
    total_repeating_non_ambiguous_double_rows = repeating_non_ambiguous_double_rows['Count'].sum()

    # Output the results for each subset
    print("\nRepeating Ambiguous Triple Rows (Non-Empty XL C):")
    print(repeating_ambiguous_triple_rows.to_markdown(index=False))
    print(f"\nTotal Repeating Ambiguous Triple Rows: {total_repeating_ambiguous_triple_rows}")

    print("\nRepeating Ambiguous Double Rows (Empty XL C):")
    print(repeating_ambiguous_double_rows.to_markdown(index=False))
    print(f"\nTotal Repeating Ambiguous Double Rows: {total_repeating_ambiguous_double_rows}")

    print("\nRepeating Non-Ambiguous Triple Rows (Non-Empty XL C):")
    print(repeating_non_ambiguous_triple_rows.to_markdown(index=False))
    print(f"\nTotal Repeating Non-Ambiguous Triple Rows: {total_repeating_non_ambiguous_triple_rows}")

    print("\nRepeating Non-Ambiguous Double Rows (Empty XL C):")
    print(repeating_non_ambiguous_double_rows.to_markdown(index=False))
    print(f"\nTotal Repeating Non-Ambiguous Double Rows: {total_repeating_non_ambiguous_double_rows}")

    # Remove character 'K' from 'XL A', 'XL B', and 'XL C' columns in the final DataFrames
    repeating_ambiguous_triple_rows['XL A'] = repeating_ambiguous_triple_rows['XL A'].str.replace('K', '')
    repeating_ambiguous_triple_rows['XL B'] = repeating_ambiguous_triple_rows['XL B'].str.replace('K', '')
    repeating_ambiguous_triple_rows['XL C'] = repeating_ambiguous_triple_rows['XL C'].str.replace('K', '')

    repeating_ambiguous_double_rows['XL A'] = repeating_ambiguous_double_rows['XL A'].str.replace('K', '')
    repeating_ambiguous_double_rows['XL B'] = repeating_ambiguous_double_rows['XL B'].str.replace('K', '')

    repeating_non_ambiguous_triple_rows['XL A'] = repeating_non_ambiguous_triple_rows['XL A'].str.replace('K', '')
    repeating_non_ambiguous_triple_rows['XL B'] = repeating_non_ambiguous_triple_rows['XL B'].str.replace('K', '')
    repeating_non_ambiguous_triple_rows['XL C'] = repeating_non_ambiguous_triple_rows['XL C'].str.replace('K', '')

    repeating_non_ambiguous_double_rows['XL A'] = repeating_non_ambiguous_double_rows['XL A'].str.replace('K', '')
    repeating_non_ambiguous_double_rows['XL B'] = repeating_non_ambiguous_double_rows['XL B'].str.replace('K', '')

    # Rename the necessary columns in the DataFrames
    repeating_ambiguous_triple_rows = repeating_ambiguous_triple_rows.rename(columns={
        'Common/Gene Name A': 'Protein1',
        'XL A': 'Residue1',
        'Common/Gene Name B': 'Protein2',
        'XL B': 'Residue2',
        'Common/Gene Name C': 'Protein3',
        'XL C': 'Residue3'
    })

    repeating_ambiguous_double_rows = repeating_ambiguous_double_rows.rename(columns={
        'Common/Gene Name A': 'Protein1',
        'XL A': 'Residue1',
        'Common/Gene Name B': 'Protein2',
        'XL B': 'Residue2'
    })

    repeating_non_ambiguous_triple_rows = repeating_non_ambiguous_triple_rows.rename(columns={
        'Common/Gene Name A': 'Protein1',
        'XL A': 'Residue1',
        'Common/Gene Name B': 'Protein2',
        'XL B': 'Residue2',
        'Common/Gene Name C': 'Protein3',
        'XL C': 'Residue3'
    })

    repeating_non_ambiguous_double_rows = repeating_non_ambiguous_double_rows.rename(columns={
        'Common/Gene Name A': 'Protein1',
        'XL A': 'Residue1',
        'Common/Gene Name B': 'Protein2',
        'XL B': 'Residue2'
    })

    header1 = ['Protein1', 'Residue1', 'Protein2', 'Residue2', 'Protein3', 'Residue3', 'Count']
    header2 = ['Protein1', 'Residue1', 'Protein2', 'Residue2', 'Count']
    # Write the final DataFrames to CSV files
    repeating_ambiguous_triple_rows.to_csv('ambiguous_triple_xls.csv', index=False, columns=header1)
    repeating_ambiguous_double_rows.to_csv('ambiguous_double_xls.csv', index=False, columns=header2)
    repeating_non_ambiguous_triple_rows.to_csv('non_ambiguous_triple_xls.csv', index=False, columns=header1)
    repeating_non_ambiguous_double_rows.to_csv('non_ambiguous_double_xls.csv', index=False, columns=header2)
    
    f1, f2 = filter_proteins(repeating_ambiguous_triple_rows, repeating_ambiguous_double_rows)
    f3, f4 = filter_proteins(repeating_non_ambiguous_triple_rows, repeating_non_ambiguous_double_rows)
    
    # Concatenate f1 with f3 and f2 with f4
    concatenated_triple = pd.concat([f1, f3])
    concatenated_double = pd.concat([f2, f4])
    print("concatenate triple", concatenated_triple)
    print("concatenate double", concatenated_double)
    
    # Remove non-unique rows
    concatenated_triple.drop(columns='Count', inplace=True)
    concatenated_double.drop(columns='Count', inplace=True)
    concatenated_triple = concatenated_triple.drop_duplicates()
    concatenated_double = concatenated_double.drop_duplicates()
    
    # Filter out rows where Protein1 and Residue1 are the same as Protein2 and Residue2
    concatenated_triple = concatenated_triple[~((concatenated_triple['Protein1'] == concatenated_triple['Protein2']) & (concatenated_triple['Residue1'] == concatenated_triple['Residue2']))]
    concatenated_double = concatenated_double[~((concatenated_double['Protein1'] == concatenated_double['Protein2']) & (concatenated_double['Residue1'] == concatenated_double['Residue2']))]
    
    # Ensure the directory exists
    output_dir = 'base_subcomplex_data'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Save the concatenated DataFrames to CSV files
    concatenated_triple.to_csv(os.path.join(output_dir, 'concatenated_triple_links.csv'), index=False)
    concatenated_double.to_csv(os.path.join(output_dir, 'concatenated_double_links.csv'), index=False)

# Replace 'full-data.csv' with your actual file path
file_path = 'full-data.csv'
analyze_repeating_rows(file_path)