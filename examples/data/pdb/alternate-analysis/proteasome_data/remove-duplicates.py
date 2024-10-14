import pandas as pd
import math
import matplotlib.pyplot as plt
from Bio.PDB import PDBParser
from scipy import stats
import json

# Function to parse PDB file and extract atom coordinates
def parse_pdb(pdb_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('structure', pdb_file)
    atom_coords = {}
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    chain_id = chain.id
                    residue_id = residue.id[1]
                    atom_name = atom.name
                    coord = atom.coord
                    atom_coords[(chain_id, residue_id, atom_name)] = coord
    return atom_coords

# Function to calculate Euclidean distance between two points
def calculate_distance(coord1, coord2):
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(coord1, coord2)))

# Function to process a DataFrame and calculate distances
def process_dataframe(df, atom_coords, output_csv, plot_title):
    distances = []
    for index, row in df.iterrows():
        chain1, residue1, atom1, chain2, residue2, atom2 = row['Chain1'], row['Residue1'], row['Atom1'], row['Chain2'], row['Residue2'], row['Atom2']
        print(chain1, residue1, atom1, chain2, residue2, atom2)
        print("atom_coords: ", atom_coords.get((chain1, residue1, atom1)))
        # try
        coord1 = atom_coords.get((chain1, residue1, atom1))
        coord2 = atom_coords.get((chain2, residue2, atom2))
        if coord1 is not None and coord2 is not None:
            distance = calculate_distance(coord1, coord2)
            distances.append(distance)
        else:
            print(f"Coordinates not found for {chain1}-{residue1}-{atom1} or {chain2}-{residue2}-{atom2}")

    # Create a DataFrame for distances
    distances_df = pd.DataFrame(distances, columns=['Distance'])

    # Save the distance distribution to a CSV file
    distances_df.to_csv(output_csv, index=False)
    print(f"Distance distribution saved to '{output_csv}'.")

    # Plot distribution of distances
    plt.figure(figsize=(10, 6))
    plt.hist(distances, bins=30, edgecolor='black', alpha=0.7)
    plt.xlabel('Distance (Å)')
    plt.ylabel('Frequency')
    plt.title(plot_title)
    plt.grid(True)
    plt.show()

    # Calculate and print statistical parameters
    mean_distance = distances_df['Distance'].mean()
    std_distance = distances_df['Distance'].std()
    max_distance = distances_df['Distance'].max()
    min_distance = distances_df['Distance'].min()
    median_distance = distances_df['Distance'].median()
    mode_distance = distances_df['Distance'].mode()[0]

    print(f"Mean Distance: {mean_distance:.2f} Å")
    print(f"Standard Deviation: {std_distance:.2f} Å")
    print(f"Maximum Distance: {max_distance:.2f} Å")
    print(f"Minimum Distance: {min_distance:.2f} Å")
    print(f"Median Distance: {median_distance:.2f} Å")
    print(f"Mode Distance: {mode_distance:.2f} Å")

    # Additional statistical parameters
    print(f"Variance: {distances_df['Distance'].var():.2f} Å^2")
    print(f"Skewness: {distances_df['Distance'].skew():.2f}")
    print(f"Kurtosis: {distances_df['Distance'].kurtosis():.2f}")
    print(f"25th Percentile: {distances_df['Distance'].quantile(0.25):.2f} Å")
    print(f"75th Percentile: {distances_df['Distance'].quantile(0.75):.2f} Å")
    print(f"Interquartile Range (IQR): {distances_df['Distance'].quantile(0.75) - distances_df['Distance'].quantile(0.25):.2f} Å")

# Read the PDB file and calculate distances
pdb_file = '5gjr.pdb'
atom_coords = parse_pdb(pdb_file)
# Convert tuple keys to strings and ndarray values to lists
atom_coords_str_keys = {str(key): value.tolist() for key, value in atom_coords.items()}

# Define the output file path
output_file = 'atom_coords.json'

# Write atom_coords to the file
with open(output_file, 'w') as file:
    json.dump(atom_coords_str_keys, file, indent=4)

print(f"atom_coords has been written to {output_file}")

# Process triples-chainid.csv
# Read the DataFrame from a CSV file
df = pd.read_csv('cleaned_triple_links.csv')

# Remove duplicate rows
df = df.drop_duplicates()

# Save the cleaned DataFrame to a new CSV file
df.to_csv('unique_triples.csv', index=False)

print("Duplicate rows removed and cleaned data saved to 'unique_triples.csv'.")
#--------------------------------------------------------------------------
# Now comes the triple links, here we will follow one code snippet from 
# another script file
#--------------------------------------------------------------------------
# Read the CSV file into a DataFrame, skipping the first row                                
headers = ['Protein1', 'Residue1', 'Protein2', 'Residue2', 'Protein3', 'Residue3']          
triple_df = pd.read_csv('unique_triples.csv', header=0, names=headers)                        
                                                                                            
# Initialize an empty list to store the new rows                                            
new_rows = []                                                                               
                                                                                            
# Process each row in the DataFrame                                                         
for index, row in triple_df.iterrows():                                                            
    new_rows.append([row['Protein1'], row['Residue1'], row['Protein2'], row['Residue2']])   
    new_rows.append([row['Protein2'], row['Residue2'], row['Protein3'], row['Residue3']])   
    new_rows.append([row['Protein3'], row['Residue3'], row['Protein1'], row['Residue1']])   
                                                                                            
# Create a new DataFrame from the new rows                                                  
new_triple_df = pd.DataFrame(new_rows, columns=['Protein1', 'Residue1', 'Protein2', 'Residue2'])
new_triple_df.to_csv('paired_triples.csv', index=False)
print("Triple links processed and saved to 'paired_triples.csv'.")
#--------------------------------------------------------------------------
# Read the csv file with the paired triples and convert it into 
# a format which biopython can read and calculate distances
#--------------------------------------------------------------------------
headers = ['Protein', 'Chain ID']
map_chainID = pd.read_csv('map-chainID.csv', names=headers)

# Replace protein names with chain IDs
new_triple_df = new_triple_df.merge(map_chainID, left_on='Protein1', right_on='Protein', how='left')
new_triple_df = new_triple_df.drop(columns=['Protein1', 'Protein'])
new_triple_df = new_triple_df.rename(columns={'Chain ID': 'Chain1'})

new_triple_df = new_triple_df.merge(map_chainID, left_on='Protein2', right_on='Protein', how='left')
new_triple_df = new_triple_df.drop(columns=['Protein2', 'Protein'])
new_triple_df = new_triple_df.rename(columns={'Chain ID': 'Chain2'})

# Add atom1 and atom2 columns with value 'CA'
new_triple_df['Atom1'] = 'CA'
new_triple_df['Atom2'] = 'CA'

# Reorder columns as chainid1, residue1, atom1, chainid2, residue2, atom2
new_triple_df = new_triple_df[['Chain1', 'Residue1', 'Atom1', 'Chain2', 'Residue2', 'Atom2']]

# Save the updated DataFrame to a new CSV file
new_triple_df.to_csv('triples-chainid.csv', index=False)
print("Protein names replaced with chain IDs and saved to 'triples-chainid.csv'.")

# Process triples DataFrame
process_dataframe(new_triple_df, atom_coords, 'exp-triple-dist.csv', 'Distribution of Calculated Distances for Triples')

# Process complete-doubles.csv
headers = ['Protein1', 'Residue1', 'Protein2', 'Residue2', 'Count']
doubles_df = pd.read_csv('complete-doubles.csv', header=0, names=headers)
doubles_df = doubles_df.drop(columns=['Count'])

# Print the DataFrame to verify the column has been dropped
print(doubles_df.head())

# Replace protein names with chain IDs for doubles
doubles_df = doubles_df.merge(map_chainID, left_on='Protein1', right_on='Protein', how='left')
doubles_df = doubles_df.drop(columns=['Protein1', 'Protein'])
doubles_df = doubles_df.rename(columns={'Chain ID': 'Chain1'})

doubles_df = doubles_df.merge(map_chainID, left_on='Protein2', right_on='Protein', how='left')
doubles_df = doubles_df.drop(columns=['Protein2', 'Protein'])
doubles_df = doubles_df.rename(columns={'Chain ID': 'Chain2'})
doubles_df = doubles_df.dropna(how='any')
print("intermediate: ", doubles_df)

# Add atom1 and atom2 columns with value 'CA'
doubles_df['Atom1'] = 'CA'
doubles_df['Atom2'] = 'CA'

# Reorder columns as chainid1, residue1, atom1, chainid2, residue2, atom2
doubles_df = doubles_df[['Chain1', 'Residue1', 'Atom1', 'Chain2', 'Residue2', 'Atom2']]
doubles_df.to_csv('doubles-chainid.csv', index=False)

# Process doubles DataFrame
process_dataframe(doubles_df, atom_coords, 'exp-double-dist.csv', 'Distribution of Calculated Distances for Doubles')
print("Triple links processed and saved to 'paired_triples.csv'.")
