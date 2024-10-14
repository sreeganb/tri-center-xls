import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import entropy
from Bio import PDB
import os
from generate_synthetic_xls import LysineCrosslinkAnalyzer  # Assuming the class is in a file named generate_synthetic_xls.py

def calculate_kl_divergence(p, q, bins=30):
    """
    Calculate the Kullback-Leibler divergence between two distributions.
    
    Parameters:
    p (array-like): First distribution.
    q (array-like): Second distribution.
    bins (int): Number of bins to use for the histogram.
    
    Returns:
    float: KL divergence.
    """
    # Estimate the distributions using histograms
    p_hist, bin_edges = np.histogram(p, bins=bins, density=True)
    q_hist, _ = np.histogram(q, bins=bin_edges, density=True)
    
    # Avoid division by zero and log of zero by adding a small constant
    p_hist += 1e-10
    q_hist += 1e-10
    
    # Normalize the histograms
    p_hist /= np.sum(p_hist)
    q_hist /= np.sum(q_hist)
    
    return entropy(p_hist, q_hist)

def check_crosslinks_between_chains(analyzer, chain1='1', chain2='y'):
    """
    Check if there are crosslinks between chainID 1 and chainID y.
    
    Parameters:
    analyzer (LysineCrosslinkAnalyzer): The analyzer object.
    chain1 (str): The first chain ID.
    chain2 (str): The second chain ID.
    
    Returns:
    bool: True if there are crosslinks between the chains, False otherwise.
    """
    crosslinks = analyzer.distances_df
    crosslinks_between_chains = crosslinks[(crosslinks['ChainID1'] == chain1) & (crosslinks['ChainID2'] == chain2)]
    return not crosslinks_between_chains.empty

def check_protein_occurrence(analyzer, protein1='Rpn2', protein2='Rpt3', min_occurrences=3):
    """
    Check if the proteins occur in the same row at least a specified number of times.
    
    Parameters:
    analyzer (LysineCrosslinkAnalyzer): The analyzer object.
    protein1 (str): The first protein name.
    protein2 (str): The second protein name.
    min_occurrences (int): Minimum number of occurrences required.
    
    Returns:
    bool: True if the proteins occur in the same row at least min_occurrences times, False otherwise.
    """
    combined_triplets = analyzer.combined_triplets
    count = combined_triplets.apply(lambda row: protein1 in row.values and protein2 in row.values, axis=1).sum()
    return count >= min_occurrences

def main():
    input_pdb = 'data/pdb/base_proteasome.pdb'
    output_pdb = 'lysine_residues.pdb'
    distance_threshold = 30  # Set the distance threshold to 30 Å
    kl_threshold = 0.3  # Set the KL divergence threshold for similarity
    reference_distribution_path = 'exp-distance-distribution.csv'  # Path to the reference distribution file

    # Load the reference distribution
    reference_distribution = pd.read_csv(reference_distribution_path)['Distance'].values

    analyzer = LysineCrosslinkAnalyzer(input_pdb, output_pdb, distance_threshold)

    while True:
        analyzer.extract_lysine_residues()
        analyzer.pick_random_crosslinks()

        # Get the generated distance distribution
        generated_distances = analyzer.distances_df['Distance'].values

        # Calculate the KL divergence
        kl_divergence = calculate_kl_divergence(generated_distances, reference_distribution)
        print(f"KL Divergence: {kl_divergence}")

        # Check if the KL divergence is below the threshold
        if kl_divergence < kl_threshold:
            print("Accepted distribution with KL divergence below threshold.")
            
            # Check for crosslinks between chainID 1 and chainID y
            if check_crosslinks_between_chains(analyzer):
                print("Crosslinks found between chainID 1 and chainID y.")
                
                # Check if Rpn2 and Rpt3 occur in the same row at least three times
                if check_protein_occurrence(analyzer):
                    print("Rpn2 and Rpt3 occur in the same row at least three times.")
                    break
                else:
                    print("Rpn2 and Rpt3 do not occur in the same row at least three times. Recalculating crosslinks...")
            else:
                print("No crosslinks found between chainID 1 and chainID y. Recalculating crosslinks...")
        else:
            print("KL divergence above threshold. Recalculating crosslinks...")

if __name__ == "__main__":
    main()