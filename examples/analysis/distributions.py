import pandas as pd
from Bio.PDB import PDBParser
import math
import matplotlib.pyplot as plt
import seaborn as sns

# Define the mapping from Protein to Chain (corrected)
# trying to be exhaustive here, get all the names possibly 
# from the experimental dataset
protein_to_chain = {
    'Rpt6': 'x',
    'Rpt3': 'y',
    'Rpt4': 'z',
    'Rpn2': '1',
    'Rpt5': '0',
    'Rpt2': 'w',
    'Rpt1': 'v',
    'Rpn11': '9',
    'Rpn3': '6',
    'Rpn5': '3',
    'Rpn9': '2',
    'alpha3': 'D',
    'alpha6': 'G',
    'alpha7': 'X',
    'beta5': 'F',
    'alpha4': 'E',
    'alpha1': 'B',
    'beta6': 'f',
    'alpha2': 'C',
    'beta4': 'd',
    'SHFM1': 'AB',
    'Rpn7': '5',
    'Rpn8': '8',
}

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

def format_data(input_file, output_file):
    # Read the input file into a DataFrame
    # data in the form of Protein1, Residue1, Protein2, Residue2
    df = pd.read_csv(input_file)

    # Use the mapping in protein_to_chain to add columns next to Protein1 and Protein2 as chain1 and chain2
    df['Chain1'] = df['Protein1'].map(protein_to_chain)
    df['Chain2'] = df['Protein2'].map(protein_to_chain)
    
    # Add two columns 'Atom1' and 'Atom2' with value 'CA'
    df['Atom1'] = 'CA'
    df['Atom2'] = 'CA'
    
    df.to_csv(output_file, index=False)

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

def plot_distance_distribution(distances, output_directory):
    # Remove None values from distances
    distances = [d for d in distances if d is not None]
    dis = pd.DataFrame(distances)
    print("mean: ", dis.mean())
    print("std: ", dis.std())
    print("len: ", len(dis))
    print("max: ", dis.max())
    plt.figure(figsize=(10, 6))
    sns.histplot(distances, bins=20)
    plt.title('Distance Distribution')
    plt.xlabel('Distance (Å)')
    plt.ylabel('Frequency')
    plt.savefig(f'{output_directory}/distance_distribution.png', dpi=300)
    plt.show()

# Main function to read CSV, parse PDB, and calculate distances
def main(pdb_file, csv_file, output_csv, is_dsso=False):
    # Parse the PDB file
    atom_coords = parse_pdb(pdb_file)
    
    formatted_csv = 'data/formatted.csv'
    
    format_data(csv_file, formatted_csv)

    # Read the CSV file into a DataFrame
    df = pd.read_csv(formatted_csv)

    # Process the DataFrame to calculate distances
    distances = process_dataframe(df, atom_coords)

    # Save the distances to a new CSV file
    distances_df = pd.DataFrame(distances, columns=['Chain1', 'Residue1', 'Atom1', 'Chain2', 'Residue2', 'Atom2', 'Distance'])
    distances_df.to_csv(output_csv, index=False)
    distances_df['Distance'].to_csv('data/exp_distances.csv', index=False)
    
    if is_dsso:
        return distances_df
    else:
        # Graphing
        # Group the distances by row number: 0,3,6,9,... is group1; 1,4,7,10,... is group2; 2,5,8,11,... is group3
        df = pd.read_csv('data/exp_distances.csv')
        df['Group'] = df.index % 3
        df1 = df[df['Group'] == 0]
        df2 = df[df['Group'] == 1]
        df3 = df[df['Group'] == 2]
        
        print("Mean of the distances for group1: ", df1['Distance'].mean())
        print("Std of the distances for group1: ", df1['Distance'].std())
        print("Mean of the distances for group2: ", df2['Distance'].mean())
        print("Std of the distances for group2: ", df2['Distance'].std())
        print("Mean of the distances for group3: ", df3['Distance'].mean())
        print("Std of the distances for group3: ", df3['Distance'].std())
        
        # Plotting the distances for df1, df2, and df3 on the same graph
        plt.figure(figsize=(10, 6))
        sns.histplot(df1['Distance'], bins=5, label='len1')
        sns.histplot(df2['Distance'], bins=5, label='len2')
        sns.histplot(df3['Distance'], bins=5, label='len3')
        
        plt.title('Distance Distribution for Triple Crosslinks')
        plt.xlabel('Distance (Å)')
        plt.ylabel('Frequency')
        plt.legend()
        plt.savefig('data/tri-distribution.png', dpi=300)
        #plt.show()
    
# Example usage
if __name__ == "__main__":
    pdb_file = 'data/5gjr.pdb'
    #csv_file = 'data/paired_tri_xls.csv'
    csv_file = 'data/dsso-xls.csv'
    #output_csv = 'data/paired_tri_distances.csv'
    output_csv = 'data/dsso_distances.csv'
    is_dsso = True
    main(pdb_file, csv_file, output_csv, is_dsso)
    if is_dsso:
        plot_distance_distribution(main(pdb_file, csv_file, output_csv, is_dsso=True)['Distance'], 'data')
