#!/usr/bin/env python3.11

import pandas as pd
import os

def ensure_directory(dir_path):
    """Ensure that a directory exists."""
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)

# Define the mapping from (Protein, Copy) to Chain
protein_copy_to_chain = {
    ('DDI1', 0): 'A',
    ('DDI1', 1): 'B',
    ('DDI2', 0): 'C',
    ('DDI2', 1): 'D',
}

ncopies = {
    'DDI1': 2,
    'DDI2': 2
}

def get_valid_copy_combinations(proteins, residues):
    """
    Determine all valid (copy1, copy2, copy3) combinations for trifunctional XL.
    Same logic as in your modeling script.
    """
    def copy_options(prot1, res1, prot2, res2):
        """Get possible copy pairs for a protein-residue pair"""
        if prot1 != prot2:
            return [(0,0), (0,1), (1,0), (1,1)]  # Different proteins: all combos
        elif res1 == res2:
            return [(0,1), (1,0)]  # Same residue: MUST be different copies
        else:
            return [(0,0), (1,1), (0,1), (1,0)]  # Different residues: ambiguous
    
    # Get possible copy pairs for each edge of the triangle
    combos_12 = copy_options(proteins[0], residues[0], proteins[1], residues[1])
    combos_23 = copy_options(proteins[1], residues[1], proteins[2], residues[2])
    combos_13 = copy_options(proteins[0], residues[0], proteins[2], residues[2])
    
    # Find consistent combinations
    valid = set()
    for (c1_12, c2_12) in combos_12:
        for (c2_23, c3_23) in combos_23:
            for (c1_13, c3_13) in combos_13:
                if c1_12 == c1_13 and c2_12 == c2_23 and c3_23 == c3_13:
                    valid.add((c1_12, c2_12, c3_23))
    
    return sorted(valid)


def process_triple_links_for_chimerax(input_file, output_file, remove_duplicates=False):
    """
    Process trifunctional crosslinks with proper ambiguity handling.
    Outputs pairwise links in ChimeraX-compatible format:
    Chain1, Residue1, Atom1, Chain2, Residue2, Atom2
    
    Parameters:
    - input_file: Path to trifunctional crosslink CSV
    - output_file: Path to output CSV
    - remove_duplicates: If True, remove duplicate pairwise links (NOT recommended 
                        for trifunctional XLs as they represent different triangular constraints)
    """
    # Read the CSV file
    headers = ['Protein1', 'Residue1', 'Protein2', 'Residue2', 'Protein3', 'Residue3']
    df = pd.read_csv(input_file, header=0, names=headers)
    
    # Initialize list to store new rows
    new_rows = []
    
    # Track statistics
    total_combos = 0
    links_per_xl = []
    
    # Process each trifunctional crosslink
    for index, row in df.iterrows():
        p1 = row['Protein1']
        p2 = row['Protein2']
        p3 = row['Protein3']
        r1 = int(''.join(filter(str.isdigit, str(row['Residue1']))))
        r2 = int(''.join(filter(str.isdigit, str(row['Residue2']))))
        r3 = int(''.join(filter(str.isdigit, str(row['Residue3']))))
        
        proteins = [p1, p2, p3]
        residues = [r1, r2, r3]
        
        # Get valid copy combinations
        valid_combos = get_valid_copy_combinations(proteins, residues)
        total_combos += len(valid_combos)
        
        print(f"\nXL {index}: {p1}:{r1} - {p2}:{r2} - {p3}:{r3}")
        print(f"  Valid copy combinations: {len(valid_combos)}")
        
        xl_links = 0
        
        # For each valid combination, add the three pairwise links
        for combo in valid_combos:
            c1, c2, c3 = combo
            
            # Get chain IDs for each copy
            chain1 = protein_copy_to_chain.get((p1, c1))
            chain2 = protein_copy_to_chain.get((p2, c2))
            chain3 = protein_copy_to_chain.get((p3, c3))
            
            if chain1 and chain2 and chain3:
                print(f"    Combo ({c1},{c2},{c3}) -> Chains ({chain1},{chain2},{chain3})")
                
                # Add three pairwise links in ChimeraX format
                # Format: Chain1, Residue1, Atom1, Chain2, Residue2, Atom2
                new_rows.append([chain1, r1, 'CA', chain2, r2, 'CA'])
                new_rows.append([chain2, r2, 'CA', chain3, r3, 'CA'])
                new_rows.append([chain3, r3, 'CA', chain1, r1, 'CA'])
                
                xl_links += 3
                
                print(f"      {chain1}:{r1} - {chain2}:{r2}")
                print(f"      {chain2}:{r2} - {chain3}:{r3}")
                print(f"      {chain3}:{r3} - {chain1}:{r1}")
        
        links_per_xl.append(xl_links)
        print(f"  Total pairwise links for this XL: {xl_links}")
    
    # Create DataFrame in ChimeraX format
    new_df = pd.DataFrame(new_rows, 
                         columns=['Chain1', 'Residue1', 'Atom1', 
                                 'Chain2', 'Residue2', 'Atom2'])
    
    print(f"\n{'='*60}")
    print(f"BEFORE DUPLICATE REMOVAL")
    print(f"{'='*60}")
    print(f"Total trifunctional XLs: {len(df)}")
    print(f"Total copy combinations: {total_combos}")
    print(f"Total pairwise links generated: {len(new_df)}")
    print(f"Expected: {total_combos * 3} (3 links per combination)")
    print(f"Links per XL: {links_per_xl}")
    
    if remove_duplicates:
        # Remove duplicates (bidirectional crosslinks are the same)
        # Sort each row to make comparisons consistent
        new_df['sorted_key'] = new_df.apply(
            lambda row: tuple(sorted([(row['Chain1'], row['Residue1'], row['Atom1']),
                                      (row['Chain2'], row['Residue2'], row['Atom2'])])),
            axis=1
        )
        initial_count = len(new_df)
        new_df = new_df.drop_duplicates(subset='sorted_key').drop(columns='sorted_key')
        removed_count = initial_count - len(new_df)
        
        print(f"\nAFTER DUPLICATE REMOVAL")
        print(f"Removed {removed_count} duplicate entries")
        print(f"Remaining unique pairwise links: {len(new_df)}")
    else:
        print(f"\nNOTE: Duplicates NOT removed - each trifunctional XL represents")
        print(f"      a unique triangular constraint even if some edges overlap")
    
    # Save to CSV (ChimeraX script expects this format)
    new_df.to_csv(output_file, index=False, header=True)
    
    print(f"\n{'='*60}")
    print(f"OUTPUT SUMMARY")
    print(f"{'='*60}")
    print(f"Processed {len(df)} trifunctional crosslinks")
    print(f"Generated {len(new_df)} pairwise links")
    print(f"Output saved to: {output_file}")
    print(f"\nFormat: Chain1, Residue1, Atom1, Chain2, Residue2, Atom2")
    print(f"Preview:\n{new_df.head(10)}")
    
    return new_df


def create_chimerax_script(crosslink_file, pdb_file, output_script='visualize_crosslinks.cxc'):
    """
    Create a ChimeraX command script to visualize crosslinks.
    
    Parameters:
    - crosslink_file: Path to the formatted crosslink CSV
    - pdb_file: Path to the PDB structure file
    - output_script: Output ChimeraX script filename
    """
    df = pd.read_csv(crosslink_file)
    
    with open(output_script, 'w') as f:
        # Header
        f.write("# ChimeraX script to visualize trifunctional crosslinks\n\n")
        
        # Load structure
        f.write(f"open {pdb_file}\n")
        f.write("cartoon\n")
        f.write("color khaki cartoons\n\n")
        
        # Color chains
        f.write("# Color chains\n")
        f.write("color /A salmon cartoons\n")
        f.write("color /B lightblue cartoons\n")
        f.write("color /C yellow cartoons\n")
        f.write("color /D green cartoons\n\n")
        
        # Show crosslinked residues
        f.write("# Show crosslinked residues as spheres\n")
        unique_residues = set()
        for _, row in df.iterrows():
            unique_residues.add((row['Chain1'], row['Residue1']))
            unique_residues.add((row['Chain2'], row['Residue2']))
        
        for chain, res in sorted(unique_residues):
            f.write(f"show /{chain}:{res}@CA atoms\n")
            f.write(f"size /{chain}:{res}@CA atomRadius 0.5\n")
        
        f.write("\n# Draw crosslinks as distance pseudobonds\n")
        
        # Draw crosslinks
        for idx, row in df.iterrows():
            chain1, res1, atom1 = row['Chain1'], row['Residue1'], row['Atom1']
            chain2, res2, atom2 = row['Chain2'], row['Residue2'], row['Atom2']
            
            # Calculate if it would be a violation (optional - shown in color)
            f.write(f"distance /{chain1}:{res1}@{atom1} /{chain2}:{res2}@{atom2}")
            f.write(f" color dark_green radius 0.6 name xl_{idx}\n")
        
        # Final settings
        f.write("\n# Final view settings\n")
        f.write("lighting soft\n")
        f.write("graphics silhouettes true\n")
        f.write("set bgColor white\n")
        f.write("view\n")
    
    print(f"\nChimeraX script created: {output_script}")
    print(f"To use in ChimeraX: open {output_script}")
    print(f"Or from command line: chimerax {output_script}")


def main():
    # Ensure output directory exists
    output_directory = 'pymol_data'
    ensure_directory(output_directory)
    
    # Process trifunctional crosslinks with ambiguity handling
    input_xl_file = 'input_data/reduced_ddi_trifunctional.csv'
    chimerax_output = os.path.join(output_directory, 'paired_double_links.csv')
    
    # Process WITHOUT removing duplicates (each triangle is unique)
    df_chimerax = process_triple_links_for_chimerax(input_xl_file, chimerax_output, 
                                                     remove_duplicates=False)
    
    # Create ChimeraX visualization script
    pdb_file = './input_data/Q8WTU0_V1_3.pdb'
    script_output = os.path.join(output_directory, 'visualize_crosslinks.cxc')
    create_chimerax_script(chimerax_output, pdb_file, script_output)
    
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    print(f"Input file: {input_xl_file}")
    print(f"ChimeraX-compatible output: {chimerax_output}")
    print(f"ChimeraX script: {script_output}")
    print("\nTo visualize in ChimeraX:")
    print(f"  Option 1 (Python): Run chimera_xls.py with visualize_crosslinks(session)")
    print(f"  Option 2 (Script): chimerax {script_output}")
    print("\nOutput format matches ChimeraX expectations:")
    print("  Chain1, Residue1, Atom1, Chain2, Residue2, Atom2")


if __name__ == "__main__":
    main()