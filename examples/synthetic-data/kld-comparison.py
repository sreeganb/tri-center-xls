import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import entropy
from Bio import PDB
import os
from generate_synthetic_xls import LysineCrosslinkAnalyzer  # Assuming the class is in a file named generate_synthetic_xls.py

def calculate_kl_divergence(p, q):
    """
    Calculate the Kullback-Leibler divergence between two distributions.
    
    Parameters:
    p (array-like): First distribution.
    q (array-like): Second distribution.
    
    Returns:
    float: KL divergence.
    """
    p = np.asarray(p, dtype=np.float)
    q = np.asarray(q, dtype=np.float)
    
    # Normalize the distributions
    p /= np.sum(p)
    q /= np.sum(q)
    
    return entropy(p, q)

def main():
    input_pdb = 'data/pdb/base_proteasome.pdb'
    output_pdb = 'lysine_residues.pdb'
    distance_threshold = 30  # Set the distance threshold to 30 Å
    kl_threshold = 0.1  # Set the KL divergence threshold for similarity
    reference_distribution_path = 'reference_distribution.csv'  # Path to the reference distribution file

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
            break
        else:
            print("KL divergence above threshold. Recalculating crosslinks...")

if __name__ == "__main__":
    main()