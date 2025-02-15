#!/usr/bin/env python3.12

import pandas as pd

def remove_k_from_columns(row):
    """
    Removes the letter 'K' from the entries in columns 2 and 4 of the row.
    
    Args:
        row (list): The input row.
        
    Returns:
        list: The modified row with 'K' removed from columns 2 and 4.
    """
    row[1] = row[1].replace('K', '')
    row[3] = row[3].replace('K', '')
    return row

def process_csv(file_path):
    """
    Processes a CSV file to remove 'K' from columns 2 and 4.
    
    Args:
        file_path (str): The path to the CSV file to be processed.
        
    Returns:
        DataFrame: A DataFrame containing the processed data.
    """
    df_raw = pd.read_csv(file_path, header=None, delimiter='\t')
    
    # Apply the removal function to each row
    processed_data = df_raw.apply(remove_k_from_columns, axis=1)
    
    return processed_data

def main():
    """
    Main function to execute the script.
    """
    # Specify the path to your CSV file
    file_path = './proteasome-data-full.csv'
    
    # Process the CSV file
    df_processed = process_csv(file_path)
    
    # Save the processed DataFrame to a new CSV file
    df_processed.to_csv('processed_proteasome_data.csv', index=False, header=False, sep='\t')
    
    # Print the resulting DataFrame
    print("Processed DataFrame:")
    print(df_processed)
    
if __name__ == "__main__":
    main()