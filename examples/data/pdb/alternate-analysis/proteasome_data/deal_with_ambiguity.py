#!/usr/bin/env python3.12

import pandas as pd

# Create a DataFrame
df = pd.read_csv('ambiguous_double_links.csv')

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
#---------------------------------------------------------------------------
# Find unique rows and their counts
unique_rows = df.value_counts()

# Convert to DataFrame for better readability
unique_rows_df = unique_rows.reset_index(name='Count')

# Calculate total number of rows in the original DataFrame
total_rows_original_df = len(df)

# Calculate total number of rows from the Count column
total_count = unique_rows_df['Count'].sum()

# Print the result
print("Unique rows with their counts:")
print(unique_rows_df)
print("\nTotal number of rows in the original DataFrame:", total_rows_original_df)
print("Total count from the Count column:", total_count)

# Check if totals match
if total_rows_original_df == total_count:
    print("The totals match.")
else:
    print("The totals do not match.")
#-----------------------------------------------------------------------------
# Drop rows that repeat and keep the first occurrence of each row onto a new DataFrame
df_unique = df.drop_duplicates(subset=['Protein1', 'Residue1', 'Protein2', 'Residue2'], keep='first')
#print(df_unique)

# Read the file named "unique_double_links.csv" into a new DataFrame
double_unique = pd.read_csv('cleaned_double_links.csv')

#---------------------------------------------------------------------------
# Find unique rows and their counts
unique_rows = double_unique.value_counts()

# Convert to DataFrame for better readability
unique_rows_df = unique_rows.reset_index(name='Count')

# Calculate total number of rows in the original DataFrame
total_rows_original_df = len(double_unique)

# Calculate total number of rows from the Count column
total_count = unique_rows_df['Count'].sum()

# Print the result
print("Unique rows with their counts:")
print(unique_rows_df)
print("\nTotal number of rows in the original DataFrame:", total_rows_original_df)
print("Total count from the Count column:", total_count)

# Check if totals match
if total_rows_original_df == total_count:
    print("The totals match.")
else:
    print("The totals do not match.")
#-----------------------------------------------------------------------------

# Calculate the number of times each unique row is repeated in the DataFrame
row_counts = df.groupby(['Protein1', 'Residue1', 'Protein2', 'Residue2']).size().reset_index(name='count')
#print("repeat counts: ", row_counts)

# Calculate the total number from the frequency of each row
total_count = row_counts['count'].sum()
#print("Total number from the frequency of each row:", total_count)

# Drop rows that are repeated in this DataFrame
double_unique = double_unique.drop_duplicates(subset=['Protein1', 'Residue1', 'Protein2', 'Residue2'], keep='first')

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
#print(unique_rows)

# combine the double_unique and unique_rows DataFrames 
combined_df = pd.concat([double_unique, unique_rows])
print("Number of rows in the combined DataFrame:", len(combined_df))

#print(combined_df)
combined_df.to_csv('combined_double_links.csv', index=False)
