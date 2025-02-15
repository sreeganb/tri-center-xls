#!/usr/bin/env python3

import pandas as pd
import os
from Bio import PDB
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Read in the CSV file
input_csv = 'data/paired_triples_with_chain_ids.csv'
df = pd.read_csv(input_csv)

# Function to calculate distance between two atoms in a PDB structure
def calculate_distance(structure, chain1, residue1, atom1, chain2, residue2, atom2):
    # Get the atoms from the residues
    atom1 = structure[0][chain1][(' ', residue1, ' ')][atom1]
    atom2 = structure[0][chain2][(' ', residue2, ' ')][atom2]

    # Calculate the distance
    distance = atom1 - atom2
    return distance

# Load the PDB structure
pdb_file = 'data/5gjr.pdb'
parser = PDB.PDBParser(QUIET=True)
structure = parser.get_structure('structure', pdb_file)

# Initialize a list to store the distances
distances = []

# Iterate over the DataFrame and calculate distances
for index, row in df.iterrows():
    residue1 = row['Residue1']
    residue2 = row['Residue2']
    chain1 = row['Chain1']
    chain2 = row['Chain2']
    atom1 = row['Atom1']
    atom2 = row['Atom2']
    
    try:
        distance = calculate_distance(structure, chain1, residue1, atom1, chain2, residue2, atom2)
        distances.append(distance)
    except KeyError as e:
        print(f"Error processing row {index}: {e}")
        distances.append(None)

# Add the distances to the DataFrame
df['Distance'] = distances

# Save the results to a new CSV file
output_csv = 'output_data/paired_triples_with_distances.csv'
df.to_csv(output_csv, index=False)

print(f"Distances calculated and saved to {output_csv}")

# Initialize lists to store the distances for each side of the triangles
side1_distances = []
side2_distances = []
side3_distances = []

# Divide the data into groups of three and order the sides by increasing distance
num_groups = 36
for i in range(num_groups):
    group = df.iloc[i*3:(i+1)*3]
    if group['Distance'].isnull().any():
        print(f"Skipping group {i+1} due to error in distance calculation.")
        continue
    group = group.sort_values(by='Distance').reset_index(drop=True)
    side1_distances.append(group.loc[0, 'Distance'])
    side2_distances.append(group.loc[1, 'Distance'])
    side3_distances.append(group.loc[2, 'Distance'])

# Calculate statistics for each side
def calculate_statistics(distances):
    mean = np.mean(distances)
    std = np.std(distances)
    max_distance = np.max(distances)
    min_distance = np.min(distances)
    count_above_30 = sum(d > 30 for d in distances)
    return mean, std, max_distance, min_distance, count_above_30

side1_stats = calculate_statistics(side1_distances)
side2_stats = calculate_statistics(side2_distances)
side3_stats = calculate_statistics(side3_distances)

# Print statistics to terminal
print(f"Side 1 - Mean: {side1_stats[0]:.2f}, Std: {side1_stats[1]:.2f}, Max: {side1_stats[2]:.2f}, Min: {side1_stats[3]:.2f}, Count > 30: {side1_stats[4]}")
print(f"Side 2 - Mean: {side2_stats[0]:.2f}, Std: {side2_stats[1]:.2f}, Max: {side2_stats[2]:.2f}, Min: {side2_stats[3]:.2f}, Count > 30: {side2_stats[4]}")
print(f"Side 3 - Mean: {side3_stats[0]:.2f}, Std: {side3_stats[1]:.2f}, Max: {side3_stats[2]:.2f}, Min: {side3_stats[3]:.2f}, Count > 30: {side3_stats[4]}")

# Write statistics to a CSV file
statistics = {
    'Side': ['Side 1', 'Side 2', 'Side 3'],
    'Mean': [side1_stats[0], side2_stats[0], side3_stats[0]],
    'Std': [side1_stats[1], side2_stats[1], side3_stats[1]],
    'Max': [side1_stats[2], side2_stats[2], side3_stats[2]],
    'Min': [side1_stats[3], side2_stats[3], side3_stats[3]],
    'Count > 30': [side1_stats[4], side2_stats[4], side3_stats[4]]
}
statistics_df = pd.DataFrame(statistics)
statistics_csv = 'output_data/statistics.csv'
statistics_df.to_csv(statistics_csv, index=False)

print(f"Statistics saved to {statistics_csv}")

# Write statistics to a LaTeX table
latex_table = statistics_df.to_latex(index=False)
latex_table_path = 'output_data/table.tex'
with open(latex_table_path, 'w') as f:
    f.write(latex_table)

print(f"LaTeX table saved to {latex_table_path}")

# Plot the distributions for each side of the triangles
plt.figure(figsize=(12, 8))

# Plot for side 1
sns.histplot(side1_distances, bins=5, color='blue', label='Side 1')

# Plot for side 2
sns.histplot(side2_distances, bins=5, color='green', label='Side 2')

# Plot for side 3
sns.histplot(side3_distances, bins=5, color='red', label='Side 3')

plt.title('Distribution of Distances for Each Side of the Triangles')
plt.xlabel('Distance (Å)')
plt.ylabel('Frequency')
plt.legend(loc='upper right')
plt.savefig('output_data/triangle_distances_distribution.png')
plt.show()