import csv
import numpy as np
import pandas as pd
from pymol import cmd
from pymol.cgo import *

#----------------------------------------------------------------------
# Generalize this code to include reading pdb files with multiple chains
# and calculating distances between atoms in different chains
# and drawing restraints between them
#----------------------------------------------------------------------

colors = {
    'green': [0.0, 1.0, 0.5],  # limegreen
    'red': [0.82, 0.0, 0.3]    # dubnium
}

radius = 0.5

def replace_protein_names(df, mapping):
    """Replaces protein names in the DataFrame with their corresponding chain IDs."""
    protein_name_columns = ['protein1', 'protein2']

    for col in protein_name_columns:
        if col in df.columns:
            df[col] = df[col].map(lambda x: mapping.get(x, x))
        else:
            print(f"Warning: Column '{col}' not found in the DataFrame.")

    # Add new columns for atom names
    df['atom1'] = 'CA'
    df['atom2'] = 'CA'

    return df

def drawRestraints(pdb_file, csv_file, selection='all', prefix='xl', threshold=30.0, atom='CA', quiet=1):
    # Load the PDB file
    cmd.load(pdb_file, 'structure')
    
    # Initialize counters
    satisfied_count = 0
    violated_count = 0
    
    # Read constraints from the CSV file
    with open(csv_file, 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            chain1 = row['chain1'].strip()
            res1 = row['Residue1'].strip()
            atom1 = row['Atom1'].strip()
            chain2 = row['chain2'].strip()
            res2 = row['Residue2'].strip()
            atom2 = row['Atom2'].strip()
            
            try:
                # Get coordinates of the atoms
                selection1 = f'{selection} and chain {chain1} and resi {res1} and name {atom1}'
                selection2 = f'{selection} and chain {chain2} and resi {res2} and name {atom2}'
                
                coords1 = cmd.get_coords(selection1)
                coords2 = cmd.get_coords(selection2)
                
                if coords1 is None:
                    print(f"Coordinates not found for atom1: {selection1}")
                if coords2 is None:
                    print(f"Coordinates not found for atom2: {selection2}")
                
                if coords1 is None or coords2 is None:
                    print(f"Skipping row due to missing coordinates: {row}")
                    continue
                
                x1, y1, z1 = coords1[0]
                x2, y2, z2 = coords2[0]
                
                d = np.linalg.norm(np.array([x2, y2, z2]) - np.array([x1, y1, z1]))
                
                if d <= float(threshold):
                    r1, g1, b1 = colors['green']
                    r2, g2, b2 = colors['green']
                    satisfied_count += 1
                else:
                    r1, g1, b1 = colors['red']
                    r2, g2, b2 = colors['red']
                    violated_count += 1
                
                cmd.load_cgo([CYLINDER, x1, y1, z1, x2, y2, z2, radius, r1, g1, b1, r2, g2, b2],
                             f'{prefix}_{res1}_{res2}_{atom}')
                
            except Exception as e:
                print(f"Error processing row {row}: {e}")
                
    cmd.group(prefix, f'{prefix}_*')
    
    # Ensure results are printed in the PyMOL console
    cmd.do(f'print "Number of restraints satisfied (<= {threshold} Å): {satisfied_count}"')
    cmd.do(f'print "Number of restraints violated (> {threshold} Å): {violated_count}"')

cmd.extend('drawRestraints', drawRestraints)
