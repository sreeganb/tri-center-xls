#!/usr/bin/env python3.12

import pandas as pd

# Create a DataFrame
df = pd.read_csv('filtered_ambiguous_double_links.csv')

# Replace ';' with '|' in the 'Residue1' and 'Residue2' columns
df['Residue1'] = df['Residue1'].str.replace(';', '|')
df['Residue2'] = df['Residue2'].str.replace(';', '|')

# Split the 'Residue1' and 'Residue2' columns by the '|' character and expand into separate rows
df['Residue1'] = df['Residue1'].str.split('|')
df['Residue2'] = df['Residue2'].str.split('|')

# Explode the DataFrame to create separate rows for each residue
df = df.explode('Residue1').explode('Residue2')

# Ensure the data types are consistent
df['Residue1'] = df['Residue1'].astype(str)
df['Residue2'] = df['Residue2'].astype(str)

# Print the resulting DataFrame
print(df)

# Drop rows that repeat and keep the first occurrence of each row onto a new DataFrame
df_unique = df.drop_duplicates(subset=['Protein1', 'Residue1', 'Protein2', 'Residue2'], keep='first')
print(df_unique)

# Read the file named "unique_double_links.csv" into a new DataFrame
double_unique = pd.read_csv('unique_double_links.csv')

# Ensure the data types are consistent in the second DataFrame
double_unique['Residue1'] = double_unique['Residue1'].astype(str)
double_unique['Residue2'] = double_unique['Residue2'].astype(str)

# Compare the rows of the two DataFrames and print the number of rows of each DataFrame as well as the number of rows that are common between the two DataFrames
common_rows = pd.merge(double_unique, df_unique, on=['Protein1', 'Residue1', 'Protein2', 'Residue2'])
print("Number of rows in double_unique:", len(double_unique))
print("Number of rows in df_unique:", len(df_unique))
print("Number of common rows:", len(common_rows))

# Select rows that are present in df_unique but not in double_unique
unique_rows = df_unique[~df_unique.isin(double_unique)].dropna()
print(unique_rows)

# combine the double_unique and unique_rows DataFrames 
combined_df = pd.concat([double_unique, unique_rows])
print(combined_df)
combined_df.to_csv('combined_double_links.csv', index=False)