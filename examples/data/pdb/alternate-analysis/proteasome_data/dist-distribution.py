import pandas as pd
from Bio.PDB import PDBParser
import math
import matplotlib.pyplot as plt
import seaborn as sns

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
    distances_df.to_csv(output_csv, index=False)
    distances_df['Distance'].to_csv('exp-bi-dist.csv', index=False)
    # graphing
    df = pd.read_csv('exp-combined.csv')
    print("mean of the distances: ", df.mean())
    print("std of the distances: ", df.std())
    
    # Assuming your data is in a DataFrame named 'df'
    sns.set_theme(style="whitegrid")  # Set the theme for a clean look

    plt.figure(figsize=(8, 5))  # Adjust figure size as needed

    # Create a histogram with density plot
    sns.histplot(data=df['Distances'], kde=False, color="steelblue", bins=30, alpha=0.7)

    # Add labels and title
    plt.xlabel('Distances', fontsize=12)
    plt.ylabel('Frequency', fontsize=12)
    plt.title('Distribution of XL Distances', fontsize=14)

    # Adjust tick font size
    plt.xticks(fontsize=11)
    plt.yticks(fontsize=11)

    # Save the figure
    plt.savefig('distance_histogram.png', dpi=300, bbox_inches='tight')
    plt.show()

# Example usage
if __name__ == "__main__":
    pdb_file = '5gjr.pdb'
    csv_file = 'doubles-chainid.csv'
    output_csv = 'doubles-distances.csv'
    main(pdb_file, csv_file, output_csv)
