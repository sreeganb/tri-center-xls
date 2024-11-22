import os
from Bio import PDB
import numpy as np
from scipy.stats import norm, skewnorm
import pandas as pd
import matplotlib.pyplot as plt

class LysineCrosslinkAnalyzer:
    def __init__(self, input_pdb, output_pdb, distance_threshold=30.0, mean=19.2, scale=9.4, probability_scaling_factor=4000, skewness=7.0, loc=6.8, skew_scale=14.0):
        self.input_pdb = input_pdb
        self.output_pdb = output_pdb
        self.distance_threshold = distance_threshold
        self.mean = mean
        self.scale = scale
        self.loc = loc
        self.skew_scale = skew_scale
        self.skewness = skewness
        self.probability_scaling_factor = probability_scaling_factor
        self.pair_distances = []  # To store precomputed pairwise distances
        # Mapping from chain IDs to protein names
        self.chain_to_protein = {
            'x': 'Rpt6',
            'y': 'Rpt3',
            'z': 'Rpt4',
            '0': 'Rpt5',
            'w': 'Rpt2',
            'v': 'Rpt1'
        }

    def normal_gaussian(self, x):
        #return norm.pdf(x, loc=self.mean, scale=self.scale)
        # Alternatively, use skewnorm if needed
        return skewnorm.pdf(x, self.skewness, loc=self.loc, scale=self.skew_scale)
    
    def extract_lysine_residues(self):
        """
        Extract lysine residues from the input PDB file.
        """
        output_directory = 'synthetic_data'
        if not os.path.exists(output_directory):
            os.makedirs(output_directory)

        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure('structure', self.input_pdb)

        io = PDB.PDBIO()
        
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

        lysine_selector = LysineSelect()
        io.set_structure(structure)
        io.save(self.output_pdb, lysine_selector)

        print(f"Number of lysine residues filtered out: {lysine_selector.lysine_count}")
        return lysine_selector.lysine_residues

    def precompute_pairwise_distances(self, lysine_residues):
        """
        Compute pairwise distances and probabilities for all unique lysine residue pairs.
        """
        self.pair_distances = []
        num_residues = len(lysine_residues)
        
        print("Precomputing pairwise distances and probabilities.")
        
        for i in range(num_residues - 1):
            ca1 = lysine_residues[i]['CA']
            for j in range(i + 1, num_residues):
                ca2 = lysine_residues[j]['CA']
                distance = ca1 - ca2
                if distance < self.distance_threshold:
                    probability = self.normal_gaussian(distance)
                    self.pair_distances.append(((i, j), distance, probability))
                    #print(f"Pair: Res {i}-{j}, Distance: {distance}, Probability: {probability}")
        
        print("Total pairs precomputed:", len(self.pair_distances))

    def select_synthetic_triplets(self, lysine_residues):
        """
        Select lysine triplets based on precomputed distances, calculating probabilities directly.
        """
        triplets = []
        num_residues = len(lysine_residues)
        
        print("Selecting triplets based on precomputed distances.")
        
        for (i, j), distance1_2, prob1_2 in self.pair_distances:
            ca1 = lysine_residues[i]['CA']
            ca2 = lysine_residues[j]['CA']
            for k in range(j + 1, num_residues):
                ca3 = lysine_residues[k]['CA']
                
                distance2_3 = ca2 - ca3
                distance3_1 = ca3 - ca1
                
                prob2_3 = self.normal_gaussian(distance2_3)
                prob3_1 = self.normal_gaussian(distance3_1)
                
                # Calculate overall selection probability
                overall_probability = prob1_2 * prob2_3 * prob3_1 * self.probability_scaling_factor
                
                # Ensure that at least two residues have different chain IDs
                chain_id1 = lysine_residues[i].get_parent().id
                chain_id2 = lysine_residues[j].get_parent().id
                chain_id3 = lysine_residues[k].get_parent().id
                
                if chain_id1 != chain_id2 or chain_id2 != chain_id3 or chain_id1 != chain_id3:
                    if overall_probability > 1e-6 and np.random.rand() < overall_probability:
                        triplet = (
                            lysine_residues[i].get_id()[1], chain_id1,
                            lysine_residues[j].get_id()[1], chain_id2,
                            lysine_residues[k].get_id()[1], chain_id3,
                            distance1_2, distance2_3, distance3_1
                        )
                        triplets.append(triplet)
                        #print("Triplet selected:", triplet)
                    
        print("Total triplets found:", len(triplets))
        
        # Create a DataFrame from triplets
        df_triplets = pd.DataFrame(triplets, columns=[
            'Residue1', 'ChainID1',
            'Residue2', 'ChainID2',
            'Residue3', 'ChainID3',
            'Distance12', 'Distance23', 'Distance31'
        ])
        
        if df_triplets.empty:
            print("No triplets found with the given parameters.")
            return df_triplets
        
        min_distance = 0.0
        max_distance = 50.0
        skewness = 7.6
        loc = 8.087
        scale = 14.797
        x_values = np.linspace(min_distance, max_distance, 1000)
        y_values = 18000*skewnorm.pdf(x_values, skewness, loc=loc, scale=scale)
        plt.plot(x_values, y_values, color='red', lw=2)
        # Write distances to CSV
        dis_df = pd.concat([df_triplets['Distance12'], df_triplets['Distance23'], df_triplets['Distance31']])
        dis_df.to_csv('synthetic_data/triplet_distances.csv', index=False)
        plt.hist(dis_df, bins=13, alpha=0.5, label='All Distances', edgecolor='black')
        plt.xlabel('Distance (Å)')
        plt.ylabel('Frequency')
        plt.title('Distribution of Distances for Selected Triplets')
        plt.legend()
        plt.show()
        
        # Map Chain IDs to Protein names
        df_triplets['Protein1'] = df_triplets['ChainID1'].map(self.chain_to_protein)
        df_triplets['Protein2'] = df_triplets['ChainID2'].map(self.chain_to_protein)
        df_triplets['Protein3'] = df_triplets['ChainID3'].map(self.chain_to_protein)
        
        # Select the required columns
        triplet = df_triplets[['Protein1', 'Residue1', 'Protein2', 'Residue2', 'Protein3', 'Residue3']]
        
        # Save the dataframe to a CSV file
        output_directory = 'synthetic_data'
        if not os.path.exists(output_directory):
            os.makedirs(output_directory)
        
        triplet.to_csv(os.path.join(output_directory, 'selected_triplets.csv'), index=False)
        
        return triplet

# Example usage
input_pdb = 'data/pdb/chopped_base_proteasome.pdb'
output_pdb = 'synthetic_data/lysine_residues.pdb'
analyzer = LysineCrosslinkAnalyzer(input_pdb, output_pdb)
lysine_residues = analyzer.extract_lysine_residues()
analyzer.precompute_pairwise_distances(lysine_residues)
triplet = analyzer.select_synthetic_triplets(lysine_residues)

# ... (previous code remains unchanged)

n = int(input("Enter the number of triplets to select: "))
min_frequency = 5  # Minimum frequency for each protein

proteins = ['Rpt1', 'Rpt2', 'Rpt3', 'Rpt4', 'Rpt5', 'Rpt6']
max_attempts = 10000  # To prevent infinite loops

if triplet.empty:
    print("No triplets available for selection.")
else:
    output_directory = f'replicates_{n}_xls'
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    for replicate_num in range(1, 6):
        selected_triplets_df = None
        for attempt in range(max_attempts):
            sampled_triplets = triplet.sample(n=n, replace=False)
            all_proteins = pd.concat([
                sampled_triplets['Protein1'],
                sampled_triplets['Protein2'],
                sampled_triplets['Protein3']
            ])
            protein_counts = all_proteins.value_counts()
            if all(protein_counts.get(protein, 0) >= min_frequency for protein in proteins):
                results = []
                pair_results = []
                distances = []
                protein_group_counter = 1
                for _, row in sampled_triplets.iterrows():
                    results.extend([
                        [row['Protein1'], row['Residue1'], f'p{protein_group_counter}', 1],
                        [row['Protein2'], row['Residue2'], f'p{protein_group_counter}', 1],
                        [row['Protein3'], row['Residue3'], f'p{protein_group_counter}', 1]
                    ])
                    pair_results.extend([
                        [row['Protein1'], row['Residue1'], row['Protein2'], row['Residue2']],
                        [row['Protein2'], row['Residue2'], row['Protein3'], row['Residue3']],
                        [row['Protein3'], row['Residue3'], row['Protein1'], row['Residue1']]
                    ])
                    protein_group_counter += 1

                pair_df = pd.DataFrame(pair_results, columns=['Protein1', 'Residue1', 'Protein2', 'Residue2'])

                residue_dict = {residue.get_id()[1]: residue for residue in lysine_residues}
                for _, pair in pair_df.iterrows():
                    res1 = residue_dict.get(pair['Residue1'])
                    res2 = residue_dict.get(pair['Residue2'])
                    if res1 is None or res2 is None:
                        continue
                    ca1 = res1['CA'].get_coord()
                    ca2 = res2['CA'].get_coord()
                    distance = np.linalg.norm(ca1 - ca2)
                    distances.append(distance)

                if distances:
                    mean_distance = np.mean(distances)
                    print(f"Replicate {replicate_num}, Attempt {attempt+1}: Mean distance = {mean_distance:.2f} Å")
                    if 16 <= mean_distance <= 21.0:
                        selected_triplets_df = sampled_triplets
                        print(f"Selected triplets for replicate {replicate_num} with mean distance {mean_distance:.2f} Å")
                        break
                    else:
                        print(f"Mean distance {mean_distance:.2f} Å not within 16-23 Å, retrying...")
                else:
                    print("No distances computed, retrying...")
            else:
                print(f"Attempt {attempt+1}: Protein frequency criteria not met, retrying...")
                continue

        if selected_triplets_df is None:
            print(f"Could not find a selection for replicate {replicate_num} after {max_attempts} attempts.")
        else:
            result_df = pd.DataFrame(results, columns=['Protein1', 'Residue1', 'Protein2', 'Residue2'])
            pair_df = pd.DataFrame(pair_results, columns=['Protein1', 'Residue1', 'Protein2', 'Residue2'])
            
            replicate_folder = os.path.join(output_directory, f'replicate_{replicate_num}')
            if not os.path.exists(replicate_folder):
                os.makedirs(replicate_folder)

            selected_triplets_df.to_csv(os.path.join(replicate_folder, 'selected_triplets.csv'), index=False)
            result_df.to_csv(os.path.join(replicate_folder, 'additional_random_lysine_triplets.csv'), index=False)
            pair_df.to_csv(os.path.join(replicate_folder, 'paired_triplets.csv'), index=False)

            def remove_two_rows_per_three(df):
                """
                For each group of three rows in the DataFrame, randomly remove two rows ensuring
                that the remaining row has 'protein1' != 'protein2'. If this condition is not met,
                repeat the process for that group until it is satisfied.
                """
                indices_to_keep = []

                for i in range(0, len(df), 3):
                    group = df.iloc[i:i+3]
                    if len(group) == 3:
                        attempts = 0
                        max_attempts = 10  # Prevent infinite loops in case of impossible conditions
                        success = False

                        while attempts < max_attempts and not success:
                            # Randomly select two indices to remove
                            remove_indices = group.sample(n=2, replace=False).index
                            # The index of the remaining row
                            keep_index = group.index.difference(remove_indices)[0]
                            remaining_row = df.loc[keep_index]

                            # Check if 'protein1' != 'protein2' in the remaining row
                            if remaining_row['Protein1'] != remaining_row['Protein2']:
                                indices_to_keep.append(keep_index)
                                success = True
                            else:
                                # Retry removal
                                attempts += 1

                        if not success:
                            # If after max_attempts the condition is not met, handle accordingly
                            # Here, we choose to skip this group or keep the row anyway
                            # You can adjust based on your requirements
                            print(f"Could not find a valid row in group starting at index {i}, proceeding without modification.")
                    else:
                        # For groups with less than three rows, keep them as is
                        indices_to_keep.extend(group.index.tolist())

                # Create the resulting DataFrame with only the kept rows
                result_df = df.loc[indices_to_keep].reset_index(drop=True)
                return result_df


            # Apply the function
            removed_df = remove_two_rows_per_three(pair_df)
            removed_df.to_csv(os.path.join(replicate_folder, 'removed_bi_xls.csv'), index=False)
