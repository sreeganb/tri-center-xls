import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import PDB
import os
from collections import defaultdict

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

    def __init__(self, input_pdb, output_pdb, distance_threshold=30, k=0.3, x0=28):
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
        plt.hist(self.distances_df['Distance'], bins=15, edgecolor='black')
        plt.title('Distribution of Distances Between C-alpha Atoms of Lysine Residues')
        plt.xlabel('Distance (Å)')
        plt.ylabel('Frequency')
        plt.savefig(os.path.join(output_directory, 'lysine_distances_distribution.png'))
        plt.show()

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

    def select_triplets_with_pairs(self, n, interacting_pairs_file):
        """
        Select n triplets at random and ensure the resulting set contains at least two crosslinks for each pair in the interacting pairs list.
        """
        output_directory = 'synthetic_data'
        if not os.path.exists(output_directory):
            os.makedirs(output_directory)

        # Read the interacting pairs from the CSV file
        interacting_pairs = pd.read_csv(interacting_pairs_file)
        interacting_pairs_set = set(tuple(x) for x in interacting_pairs.values)

        while True:
            selected_triplets = self.triplets_df.sample(n=n, random_state=np.random.randint(0, 10000))
            pair_counts = defaultdict(int)

            for _, row in selected_triplets.iterrows():
                pairs = [
                    (row['Residue1'], row['Residue2']),
                    (row['Residue2'], row['Residue3']),
                    (row['Residue3'], row['Residue1'])
                ]
                for pair in pairs:
                    if pair in interacting_pairs_set or (pair[1], pair[0]) in interacting_pairs_set:
                        pair_counts[pair] += 1

            # Check if each pair in the interacting pairs list appears at least two times
            valid = all(count >= 12 for pair, count in pair_counts.items() if pair in interacting_pairs_set or (pair[1], pair[0]) in interacting_pairs_set)
            if valid:
                break

        selected_triplets = selected_triplets.drop(columns=['Distance12', 'Distance23', 'Distance31'])
        selected_triplets['Protein1'] = selected_triplets['ChainID1'].map({v: k for k, v in self.protein_to_chain.items()})
        selected_triplets['Protein2'] = selected_triplets['ChainID2'].map({v: k for k, v in self.protein_to_chain.items()})
        selected_triplets['Protein3'] = selected_triplets['ChainID3'].map({v: k for k, v in self.protein_to_chain.items()})
        selected_triplets = selected_triplets[['Protein1', 'Residue1', 'Protein2', 'Residue2', 'Protein3', 'Residue3']]
        
        # Save the selected triplets to a file
        selected_triplets.to_csv(os.path.join(output_directory, 'selected_lysine_triplets.csv'), index=False)
        print(f"Selected {n} triplets with at least two crosslinks for each pair in the interacting pairs list.")
        
        # Process DataFrame in chunks of 1 row
        results = []
        pair_results = []
        protein_group_counter = 1
        for start in range(0, len(selected_triplets), 1):
            chunk = selected_triplets.iloc[start:start+1]

            # Collect unique (protein, res) pairs and their occurrences
            for _, row in chunk.iterrows():
                results.append([row['Protein1'], row['Residue1'], f'p{protein_group_counter}', 1])
                results.append([row['Protein2'], row['Residue2'], f'p{protein_group_counter}', 1])
                results.append([row['Protein3'], row['Residue3'], f'p{protein_group_counter}', 1])
               
                pair_results.append([row['Protein1'], row['Residue1'], row['Protein2'], row['Residue2']])
                pair_results.append([row['Protein2'], row['Residue2'], row['Protein3'], row['Residue3']])
                pair_results.append([row['Protein3'], row['Residue3'], row['Protein1'], row['Residue1']])

            protein_group_counter += 1

        # Create DataFrame from results
        result_df = pd.DataFrame(results, columns=['Protein1', 'Residue1', 'Protein2', 'Residue2'])
        result_df.to_csv(os.path.join(output_directory, 'additional_random_lysine_triplets.csv'), index=False)
        result_df = pd.DataFrame(pair_results, columns=['Protein1', 'Residue1', 'Protein2', 'Residue2'])
        result_df.to_csv(os.path.join(output_directory, 'paired_triplets.csv'), index=False)

# Example usage
input_pdb = 'data/pdb/base_proteasome.pdb'
output_pdb = 'lysine_residues.pdb'
distance_threshold = 30  # Set the distance threshold to 30 Å
analyzer = LysineCrosslinkAnalyzer(input_pdb, output_pdb, distance_threshold)
analyzer.extract_lysine_residues()
analyzer.select_triplets_with_pairs(n=100, interacting_pairs_file='input_data/interacting_pairs.csv')