#!/usr/bin/env python3.12

import pandas as pd
import matplotlib.pyplot as plt
import os

def replace_and_split_columns(df, columns, old, new):
    df[columns] = df[columns].apply(lambda x: x.str.replace(old, new).str.split(new))
    return df

def explode_dataframe(df, columns):
    for column in columns:
        df = df.explode(column)
    return df

def ensure_data_types(df, columns, dtype):
    df[columns] = df[columns].astype(dtype)
    return df

def find_unique_rows_and_counts(df):
    unique_rows = df.value_counts().reset_index(name='Count')
    total_rows = len(df)
    total_count = unique_rows['Count'].sum()
    return unique_rows, total_rows, total_count

def drop_duplicates(df, subset):
    return df.drop_duplicates(subset=subset, keep='first')

# Read the CSV files into DataFrames
ambiguous_df = pd.read_csv('ambiguous_double_links.csv')
cleaned_df = pd.read_csv('cleaned_double_links.csv')

# Replace ';' with '|' and split the columns in one step
ambiguous_df = replace_and_split_columns(ambiguous_df, ['Residue1', 'Residue2'], ';', '|')

# Explode the DataFrame to create separate rows for each residue
ambiguous_df = explode_dataframe(ambiguous_df, ['Residue1', 'Residue2'])

# Ensure the data types are consistent
ambiguous_df = ensure_data_types(ambiguous_df, ['Residue1', 'Residue2'], str)
cleaned_df = ensure_data_types(cleaned_df, ['Residue1', 'Residue2'], str)

# Find unique rows and their counts in ambiguous_df
unique_rows_ambiguous, total_rows_ambiguous, total_count_ambiguous = find_unique_rows_and_counts(ambiguous_df)

# Print the result for ambiguous_df
print("Unique rows with their counts in ambiguous_df:")
print(unique_rows_ambiguous)
print("\nTotal number of rows in the original ambiguous_df:", total_rows_ambiguous)
print("Total count from the Count column:", total_count_ambiguous)

# Check if totals match for ambiguous_df
if total_rows_ambiguous == total_count_ambiguous:
    print("The totals match for ambiguous_df.")
else:
    print("The totals do not match for ambiguous_df.")

# Drop duplicate rows and keep the first occurrence in ambiguous_df
unique_ambiguous_df = drop_duplicates(ambiguous_df, ['Protein1', 'Residue1', 'Protein2', 'Residue2'])

# Find unique rows and their counts in cleaned_df
unique_rows_cleaned, total_rows_cleaned, total_count_cleaned = find_unique_rows_and_counts(cleaned_df)

# Print the result for cleaned_df
print("Unique rows with their counts in cleaned_df:")
print(unique_rows_cleaned)
print("\nTotal number of rows in the original cleaned_df:", total_rows_cleaned)
print("Total count from the Count column:", total_count_cleaned)

# Check if totals match for cleaned_df
if total_rows_cleaned == total_count_cleaned:
    print("The totals match for cleaned_df.")
else:
    print("The totals do not match for cleaned_df.")

# Drop duplicate rows and keep the first occurrence in cleaned_df
unique_cleaned_df = drop_duplicates(cleaned_df, ['Protein1', 'Residue1', 'Protein2', 'Residue2'])

#*******************************************************************************
# Combine the unique ambiguous and cleaned DataFrames into a single dataframe
# and print the number of rows in the combined DataFrame
#*******************************************************************************
combined_df = pd.concat([unique_ambiguous_df, unique_cleaned_df], ignore_index=True)
print("\nNumber of rows in the combined DataFrame:", len(combined_df))
# find out the unique non repeating rows in combined_df
unique_combined, total_rows_combined, total_count_combined = find_unique_rows_and_counts(combined_df)
print("Unique rows with their counts in combined_df:")
print(unique_combined)
print("\nTotal number of rows in the original combined_df:", total_rows_combined)
print("Total count from the Count column:", total_count_combined)

# Check if totals match for combined_df
if total_rows_combined == total_count_combined:
    print("The totals match for combined_df.")
else:
    print("The totals do not match for combined_df.")

unique_combined_df = drop_duplicates(combined_df, ['Protein1', 'Residue1', 'Protein2', 'Residue2'])

# write the combined_df to a csv file
unique_combined_df.to_csv('combined_double_links.csv', index=False)

#--------------------------------------------------------------------------
# Now comes the triple links, here we will follow one code snippet from 
# another script file
#--------------------------------------------------------------------------
# Read the CSV file into a DataFrame, skipping the first row                                
headers = ['Protein1', 'Residue1', 'Protein2', 'Residue2', 'Protein3', 'Residue3']          
triple_df = pd.read_csv('triple_links.csv', header=0, names=headers)                        
                                                                                            
# Initialize an empty list to store the new rows                                            
new_rows = []                                                                               
                                                                                            
# Process each row in the DataFrame                                                         
for index, row in triple_df.iterrows():                                                            
    new_rows.append([row['Protein1'], row['Residue1'], row['Protein2'], row['Residue2']])   
    new_rows.append([row['Protein2'], row['Residue2'], row['Protein3'], row['Residue3']])   
    new_rows.append([row['Protein3'], row['Residue3'], row['Protein1'], row['Residue1']])   
                                                                                            
# Create a new DataFrame from the new rows                                                  
new_triple_df = pd.DataFrame(new_rows, columns=['Protein1', 'Residue1', 'Protein2', 'Residue2'])
# find out unique rows and their counts
unique_triple, total_rows_triple, total_count_triple = find_unique_rows_and_counts(new_triple_df)
print("Unique rows with their counts in new_triple_df:")
print(unique_triple)
print("\nTotal number of rows in the original new_triple_df:", total_rows_triple)
print("Total count from the Count column:", total_count_triple)

# Check if totals match for new_triple_df
if total_rows_triple == total_count_triple:
    print("The totals match for new_triple_df.")
else:
    print("The totals do not match for new_triple_df.")

unique_triple_df = drop_duplicates(new_triple_df, ['Protein1', 'Residue1', 'Protein2', 'Residue2'])
   
unique_triple_df.to_csv('paired_triple_links.csv', index=False)                                       
print("triple links: ", unique_triple_df)
#--------------------------------------------------------------------------
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
                                                                                                                         
    #filtered_triple_links = df_triple_links[                                                                             
    #    df_triple_links[['Protein1', 'Protein2', 'Protein3']].apply(lambda x: set(x).issubset(allowed_proteins), axis=1) 
    #]                                                                                                                    
                                                                                                                         
    filtered_double_links = df_double_links[                                                                             
        df_double_links[['Protein1', 'Protein2']].apply(lambda x: set(x).issubset(allowed_proteins), axis=1)             
    ]                                                                                                                    
                                                                                                                         
    # Ensure the directory exists                                                                                        
    output_dir = 'base_subcomplex_data'                                                                                  
    if not os.path.exists(output_dir):                                                                                   
        os.makedirs(output_dir)                                                                                          
                                                                                                                         
    # Save filtered DataFrames to CSV files    
    #print("filtered triple links: ", filtered_triple_links)
    print("filtered double links: ", filtered_double_links)                                                                          
    #filtered_triple_links.to_csv(os.path.join(output_dir, 'filtered-triple-links.csv'), index=False)                     
    filtered_double_links.to_csv(os.path.join(output_dir, 'filtered-double-links.csv'), index=False)                     
                                                                                                                         
    #return filtered_triple_links, filtered_double_links                                                                  
#--------------------------------------------------------------------------
filter_proteins(unique_triple_df, unique_combined_df)
                                                                       