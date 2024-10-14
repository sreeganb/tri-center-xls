import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import entropy
from Bio import PDB
import os
from generate_synthetic_xls import LysineCrosslinkAnalyzer  # Assuming the class is in a file named generate_synthetic_xls.py

def calculate_kl_divergence(selected_distances, reference_distribution, bins=30):
    """
    Calculate the Kullback-Leibler divergence between two distributions.
    
    Parameters:
    selected_distances (array-like): Generated distance distribution.
    reference_distribution (array-like): Reference distribution.
    bins (int): Number of bins to use for the histogram.
    
    Returns:
    float: KL divergence.
    """
    # Estimate the distributions using histograms
    p_hist, bin_edges = np.histogram(selected_distances, bins=bins, density=True)
    q_hist, _ = np.histogram(reference_distribution, bins=bin_edges, density=True)
    
    # Avoid division by zero and log of zero by adding a small constant
    p_hist += 1e-10
    q_hist += 1e-10
    
    # Normalize the histograms
    p_hist /= np.sum(p_hist)
    q_hist /= np.sum(q_hist)
    
    return entropy(p_hist, q_hist)

def check_protein_occurrence(combined_triplets, protein1='Rpn2', protein2='Rpt3', min_occurrences=3):
    """
    Check if the proteins occur in the same row at least a specified number of times.
    
    Parameters:
    combined_triplets (DataFrame): DataFrame containing the combined triplets.
    protein1 (str): The first protein name.
    protein2 (str): The second protein name.
    min_occurrences (int): Minimum number of occurrences required.
    
    Returns:
    bool: True if the proteins occur in the same row at least min_occurrences times, False otherwise.
    """
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

        # Get the generated distance distribution from selected_distances
        selected_distances = analyzer.selected_distances

        # Calculate the KL divergence
        kl_divergence = calculate_kl_divergence(selected_distances, reference_distribution)
        print(f"KL Divergence: {kl_divergence}")

        # Check if the KL divergence is below the threshold
        if kl_divergence < kl_threshold:
            print("Accepted distribution with KL divergence below threshold.")
            
            # Check if Rpn2 and Rpt3 occur in the same row at least three times
            if check_protein_occurrence(analyzer.combined_triplets):
                print("Rpn2 and Rpt3 occur in the same row at least three times.")
                break
            else:
                print("Rpn2 and Rpt3 do not occur in the same row at least three times. Recalculating crosslinks...")
        else:
            print("KL divergence above threshold. Recalculating crosslinks...")

if __name__ == "__main__":
    main()