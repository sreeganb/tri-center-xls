import pandas as pd
from Bio.PDB import PDBParser
import math
import matplotlib.pyplot as plt
from scipy.stats import skewnorm
import numpy as np

# Function to parse PDB file and extract atom coordinates
def parse_pdb(pdb_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('structure', pdb_file)
    atom_coords = {}
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    chain_id = chain.id.strip()
                    residue_id = residue.id[1]
                    atom_name = atom.name.strip()
                    coord = atom.coord
                    atom_coords[(chain_id, residue_id, atom_name)] = coord
    return atom_coords

# Function to calculate Euclidean distance between two points
def calculate_distance(coord1, coord2):
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(coord1, coord2)))

# Function to process a DataFrame and calculate distances
def process_dataframe(df, atom_coords):
    distances = []
    for index, row in df.iterrows():
        chain1 = str(row['Chain1']).strip()
        residue1 = int(row['Residue1'])
        atom1 = str(row['Atom1']).strip()
        chain2 = str(row['Chain2']).strip()
        residue2 = int(row['Residue2'])
        atom2 = str(row['Atom2']).strip()
        
        coord1 = atom_coords.get((chain1, residue1, atom1))
        coord2 = atom_coords.get((chain2, residue2, atom2))
        if coord1 is not None and coord2 is not None:
            distance = calculate_distance(coord1, coord2)
            distances.append((chain1, residue1, atom1, chain2, residue2, atom2, distance))
        else:
            print(f"Coordinates not found for: {chain1}-{residue1}-{atom1} or {chain2}-{residue2}-{atom2}")
    return distances

# Main function to read CSV, parse PDB, and calculate distances
def main(pdb_file, csv_file, output_csv):
    # Parse the PDB file
    atom_coords = parse_pdb(pdb_file)

    # Read the CSV file into a DataFrame
    df = pd.read_csv(csv_file)

    # Process the DataFrame to calculate distances
    distances = process_dataframe(df, atom_coords)

    # Save the distances to a new CSV file
    distances_df = pd.DataFrame(distances, columns=['Chain1', 'Residue1', 'Atom1', 'Chain2', 'Residue2', 'Atom2', 'Distance'])
    distances_df.drop_duplicates(inplace=True)
    distances_df.to_csv(output_csv, index=False)
    distances_df['Distance'].to_csv('exp-dist.csv', index=False)

    # Filter out distances greater than 100.0
    filtered_distances_df = distances_df[distances_df['Distance'] <= 100.0]

    mean_distance = filtered_distances_df['Distance'].mean()
    std_distance = filtered_distances_df['Distance'].std()
    max_distance = filtered_distances_df['Distance'].max()
    min_distance = filtered_distances_df['Distance'].min()

    print("mean of the distances: ", mean_distance)
    print("std of the distances: ", std_distance)
    print("max of the distances: ", max_distance)
    print("min of the distances: ", min_distance)
    
    # Fit a skewed Gaussian distribution to the data
    skewness, loc, scale = skewnorm.fit(filtered_distances_df['Distance'])
    print(f"Fitted skewness: {skewness}, loc: {loc}, scale: {scale}")
    x_values = np.linspace(min_distance, max_distance, 1000)
    y_values = skewnorm.pdf(x_values, skewness, loc=loc, scale=scale)

    # graphing
    plt.figure(figsize=(8, 5))  # Adjust figure size as needed

    # Create a histogram
    plt.hist(filtered_distances_df['Distance'], bins=21, color="steelblue", alpha=0.6, edgecolor='black', density=True)

    # Plot the skewed Gaussian distribution
    plt.plot(x_values, y_values, color='red', lw=2)

    # Add labels and title
    plt.xlabel('Distances', fontsize=12)
    plt.ylabel('Density', fontsize=12)
    plt.title('Distribution of TSTO XL Distances', fontsize=14)

    # Adjust tick font size
    plt.xticks(fontsize=11)
    plt.yticks(fontsize=11)

    # Save the figure
    plt.savefig('distance_histogram.png', dpi=300, bbox_inches='tight')
    plt.show()

# Example usage
if __name__ == "__main__":
    pdb_file = '5gjr.pdb'
    csv_file = 'xls-chainid.csv'
    output_csv = 'xls-distances.csv'
    main(pdb_file, csv_file, output_csv)