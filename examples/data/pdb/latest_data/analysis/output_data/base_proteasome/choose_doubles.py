import os
import pandas as pd
import numpy as np

def main():
    # Read the CSV files
    combined_double_links = pd.read_csv('combined_double_links.csv')
    paired_triple_links = pd.read_csv('paired_triple_links.csv')
    
    # Use 'Protein1' and 'Protein2' as subunit identifiers
    subunit_col1 = 'Protein1'
    subunit_col2 = 'Protein2'
    
    # Separate inter and intra subunit linkages
    inter_subunit_links = combined_double_links[
        combined_double_links[subunit_col1] != combined_double_links[subunit_col2]
    ]
    
    intra_subunit_links = combined_double_links[
        combined_double_links[subunit_col1] == combined_double_links[subunit_col2]
    ]
    
    n_samples = 83 - 18  # Number of rows to select
    n_reps = 5  # Number of repetitions
    min_inter_subunit = int(np.ceil(0.95 * n_samples))  # Minimum inter subunit links
    max_intra_subunit = n_samples - min_inter_subunit   # Remaining intra subunit links
    
    output_dir = 'resampling'
    os.makedirs(output_dir, exist_ok=True)
    
    for i in range(1, n_reps + 1):
        # Sample inter subunit links
        sampled_inter = inter_subunit_links.sample(
            n=min_inter_subunit, replace=False, random_state=i
        )
        
        # Sample intra subunit links
        sampled_intra = intra_subunit_links.sample(
            n=max_intra_subunit, replace=False, random_state=100 + i  # Different seed
        )
        
        # Combine sampled inter and intra subunit links
        sampled_doubles = pd.concat([sampled_inter, sampled_intra], ignore_index=True)
        
        # Shuffle the combined sample
        sampled_doubles = sampled_doubles.sample(
            frac=1, random_state=200 + i
        ).reset_index(drop=True)
        
        # Save the sampled doubles to a CSV file
        output_filename_doubles = f'rep{i}_doubles.csv'
        sampled_doubles.to_csv(os.path.join(output_dir, output_filename_doubles), index=False)
        
        # Save the triples to a separate CSV file (unchanged)
        output_filename_triples = f'rep{i}_triples.csv'
        paired_triple_links.to_csv(os.path.join(output_dir, output_filename_triples), index=False)
        
        print(f"Created {output_filename_doubles} with {len(sampled_doubles)} rows.")
        print(f"Created {output_filename_triples} with {len(paired_triple_links)} rows.")
        print(f"Inter subunit linkages in doubles: {len(sampled_inter)}")
        print(f"Intra subunit linkages in doubles: {len(sampled_intra)}\n")

if __name__ == '__main__':
    main()
