#!/usr/bin/env python3.11

import pandas as pd

# Read the CSV file into a DataFrame, skipping the first row
headers = ['Protein1', 'Residue1', 'Protein2', 'Residue2', 'Protein3', 'Residue3']
df = pd.read_csv('unique_triple_links.csv', header=0, names=headers)

# Initialize an empty list to store the new rows
new_rows = []

# Process each row in the DataFrame
for index, row in df.iterrows():
    new_rows.append([row['Protein1'], row['Residue1'], row['Protein2'], row['Residue2']])
    new_rows.append([row['Protein2'], row['Residue2'], row['Protein3'], row['Residue3']])
    new_rows.append([row['Protein3'], row['Residue3'], row['Protein1'], row['Residue1']])

# Create a new DataFrame from the new rows
new_df = pd.DataFrame(new_rows, columns=['Protein1', 'Residue1', 'Protein2', 'Residue2'])
new_df.to_csv('paired_double_links.csv', index=False)

# Read the combined_double_links.csv file into a DataFrame
combined_df = pd.read_csv('combined_double_links.csv')

# Initialize a list to store rows from df1 that match any row in df2
matching_rows = []

# Iterate over each row in df1
for index1, row1 in new_df.iterrows():
    # Check if the row exists in df2
    if row1.tolist() in combined_df.values.tolist():
        matching_rows.append(row1)

# Convert the list of matching rows to a DataFrame
df3 = pd.DataFrame(matching_rows)

# Save the DataFrame to a CSV file
#df3.to_csv('matching_rows.csv', index=False)

#print("Matching rows saved to 'matching_rows.csv'")
print(df3)