import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import PDB
import os

class LysineCrosslinkAnalyzer:
    # Define the protein to chain mapping
    protein_to_chain = {
        'Rpt6': 'x',
        'Rpt3': 'y',
        'Rpt4': 'z',
        'Rpn2': '1',
        'Rpt5': '0',
        'Rpt2': 'w',
        'Rpt1': 'v'
    }

    def __init__(self, input_pdb, output_pdb, distance_threshold=30, k=0.4, x0=28):
        """
        Initialize the LysineCrosslinkAnalyzer with input PDB file, output PDB file, and distance threshold.
        """
        self.input_pdb = input_pdb
        self.output_pdb = output_pdb
        self.distance_threshold = distance_threshold
        self.k = k
        self.x0 = x0
        self.distances_df = None
        self.triplets_df = None
        self.selected_distances = []  # Add this line
        self.combined_triplets = None  # Add this line

    def logistic_function(self, x):
        """
        Logistic function that transitions from 1 to 0.
        
        Parameters:
        x (array-like): Input values.
        
        Returns:
        array-like: Output values of the logistic function.
        """
        return 1 / (1 + np.exp(self.k * (x - self.x0)))

    def extract_lysine_residues(self):
        """
        Extract lysine residues from the input PDB file, calculate distances between them, and find triplets.
        """
        output_directory = 'synthetic_data'
        # Ensure the output directory exists
        if not os.path.exists(output_directory):
            os.makedirs(output_directory)

        # Create a PDB parser
        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure('structure', self.input_pdb)

        # Create a PDB writer
        io = PDB.PDBIO()
        
        # Define a class to select only lysine residues and count them
        class LysineSelect(PDB.Select):
            def __init__(self):
                self.lysine_count = 0
                self.lysine_residues = []

            def accept_residue(self, residue):
                if residue.get_resname() == 'LYS':
                    self.lysine_count += 1
                    self.lysine_residues.append(residue)
                    return True
                return False

        # Create an instance of LysineSelect
        lysine_selector = LysineSelect()

        # Save the lysine residues to a new PDB file
        io.set_structure(structure)
        io.save(self.output_pdb, lysine_selector)

        # Print the number of lysine residues
        print(f"Number of lysine residues filtered out: {lysine_selector.lysine_count}")

        # Calculate distances between C-alpha atoms of all pairs of lysine residues
        distances = []
        for i, res1 in enumerate(lysine_selector.lysine_residues):
            for j, res2 in enumerate(lysine_selector.lysine_residues):
                if i < j:
                    ca1 = res1['CA']
                    ca2 = res2['CA']
                    distance = ca1 - ca2
                    probability = self.logistic_function(distance)
                    if np.random.rand() < probability:  # Select based on logistic probability
                        distances.append((res1.get_id()[1], res1.get_parent().id, res2.get_id()[1], res2.get_parent().id, distance))

        # Save distances to a file
        self.distances_df = pd.DataFrame(distances, columns=['Residue1', 'ChainID1', 'Residue2', 'ChainID2', 'Distance'])
        self.distances_df.to_csv(os.path.join(output_directory, 'lysine_distances.csv'), index=False)

        # Plot the distance distribution
        plt.hist(self.distances_df['Distance'], bins=30, edgecolor='black')
        plt.title('Distribution of Distances Between C-alpha Atoms of Lysine Residues')
        plt.xlabel('Distance (Å)')
        plt.ylabel('Frequency')
        plt.savefig(os.path.join(output_directory, 'lysine_distances_distribution.png'))
        #plt.show()

        # Find triplets of lysine residues
        triplets = []
        for i, res1 in enumerate(lysine_selector.lysine_residues):
            for j, res2 in enumerate(lysine_selector.lysine_residues):
                if i < j:
                    ca1 = res1['CA']
                    ca2 = res2['CA']
                    distance1_2 = ca1 - ca2
                    probability1_2 = self.logistic_function(distance1_2)
                    if np.random.rand() < probability1_2:  # Select based on logistic probability
                        for k, res3 in enumerate(lysine_selector.lysine_residues):
                            if j < k:
                                ca3 = res3['CA']
                                distance2_3 = ca2 - ca3
                                distance3_1 = ca3 - ca1
                                probability2_3 = self.logistic_function(distance2_3)
                                probability3_1 = self.logistic_function(distance3_1)
                                if np.random.rand() < probability2_3 and np.random.rand() < probability3_1:
                                    triplets.append((res1.get_id()[1], res1.get_parent().id, res2.get_id()[1], res2.get_parent().id, res3.get_id()[1], res3.get_parent().id, distance1_2, distance2_3, distance3_1))

        # Save triplets to a file
        self.triplets_df = pd.DataFrame(triplets, columns=['Residue1', 'ChainID1', 'Residue2', 'ChainID2', 'Residue3', 'ChainID3', 'Distance12', 'Distance23', 'Distance31'])
        self.triplets_df.to_csv(os.path.join(output_directory, 'lysine_triplets.csv'), index=False)

        # Print the number of triplets found
        print(f"Number of triplets found: {len(triplets)}")

    def pick_random_crosslinks(self):
        """
        Randomly pick 30 crosslinks from the doubles and triplets, process them, and save the results to files.
        """
        # Ensure the output directory exists
        output_directory = 'synthetic_data'
        if not os.path.exists(output_directory):
            os.makedirs(output_directory)
        self.selected_distances = []  # Initialize the attribute
            
        # Randomly pick 30 crosslinks from the doubles
        random_doubles = self.distances_df.sample(n=40, random_state=9182791)
        self.selected_distances = random_doubles['Distance'].tolist()
        random_doubles = random_doubles.drop(columns=['Distance'])
        random_doubles['Protein1'] = random_doubles['ChainID1'].map({v: k for k, v in self.protein_to_chain.items()})
        random_doubles['Protein2'] = random_doubles['ChainID2'].map({v: k for k, v in self.protein_to_chain.items()})
        random_doubles = random_doubles[['Protein1', 'Residue1', 'Protein2', 'Residue2']]
        random_doubles.to_csv(os.path.join(output_directory, 'random_lysine_doubles.csv'), index=False)

        # Load the experimental triplets from the CSV file
        exp_triplets_path = os.path.join(output_directory, 'exp-triple-xls.csv')
        exp_triplets = pd.read_csv(exp_triplets_path)

        # Randomly pick the remaining crosslinks from the triplets
        remaining_triplets_count = 30 - len(exp_triplets)
        random_triplets = self.triplets_df.sample(n=remaining_triplets_count, random_state=2129081000)
        self.selected_distances.extend(random_triplets[['Distance12', 'Distance23', 'Distance31']].values.flatten().tolist())
        random_triplets = random_triplets.drop(columns=['Distance12', 'Distance23', 'Distance31'])
        random_triplets['Protein1'] = random_triplets['ChainID1'].map({v: k for k, v in self.protein_to_chain.items()})
        random_triplets['Protein2'] = random_triplets['ChainID2'].map({v: k for k, v in self.protein_to_chain.items()})
        random_triplets['Protein3'] = random_triplets['ChainID3'].map({v: k for k, v in self.protein_to_chain.items()})
        random_triplets = random_triplets[['Protein1', 'Residue1', 'Protein2', 'Residue2', 'Protein3', 'Residue3']]

        # Combine the experimental triplets with the randomly picked triplets
        self.combined_triplets = pd.concat([exp_triplets, random_triplets], ignore_index=True)  # Add this line
        self.combined_triplets.to_csv(os.path.join(output_directory, 'random_lysine_triplets.csv'), index=False)

        # Process DataFrame in chunks of 1 row
        results = []
        protein_group_counter = 1
        for start in range(0, len(self.combined_triplets), 1):
            chunk = self.combined_triplets.iloc[start:start+1]

            # Collect unique (protein, res) pairs and their occurrences
            for _, row in chunk.iterrows():
                results.append([row['Protein1'], row['Residue1'], f'p{protein_group_counter}', 1])
                results.append([row['Protein2'], row['Residue2'], f'p{protein_group_counter}', 1])
                results.append([row['Protein3'], row['Residue3'], f'p{protein_group_counter}', 1])

            protein_group_counter += 1

        # Create DataFrame from results
        result_df = pd.DataFrame(results, columns=['Protein1', 'Residue1', 'Protein2', 'Residue2'])
        result_df.to_csv(os.path.join(output_directory, 'additional_random_lysine_triplets.csv'), index=False)

        print("Randomly picked 30 crosslinks from doubles and 30 from triplets and saved to files in synthetic_data.")
        print("Generated additional file with specific format and saved to synthetic_data.")

        # Plot the distance distribution of selected pairs and triplets
        print(f"Number of selected distances: {len(self.selected_distances)}")
        print(f"mean distance: {np.mean(self.selected_distances):.4f} Å")
        print(f"std distance: {np.std(self.selected_distances):.4f} Å")
        plt.hist(self.selected_distances, bins=30, edgecolor='black')
        plt.title('Distribution of Selected Distances Between C-alpha Atoms of Lysine Residues')
        plt.xlabel('Distance (Å)')
        plt.ylabel('Frequency')
        plt.savefig(os.path.join(output_directory, 'selected_lysine_distances_distribution.png'))
        #plt.show()

# Example usage
input_pdb = 'data/pdb/base_proteasome.pdb'
output_pdb = 'lysine_residues.pdb'
distance_threshold = 30  # Set the distance threshold to 30 Å
analyzer = LysineCrosslinkAnalyzer(input_pdb, output_pdb, distance_threshold)
analyzer.extract_lysine_residues()
analyzer.pick_random_crosslinks()