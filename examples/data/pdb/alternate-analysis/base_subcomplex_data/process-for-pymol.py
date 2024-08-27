#!/usr/bin/env python3.12

import csv
import pandas as pd

# Define the mapping from Protein to Chain (corrected)
protein_to_chain = {
    'Rpt6': 'x',
    'Rpt3': 'y',
    'Rpt4': 'z',
    'Rpn2': '1',
    'Rpt5': '0',
    'Rpt2': 'w',
    'Rpt1': 'v'
}

# Input and output file names
#input_file = 'paired_triple_links.csv'
#output_file = 'converted_links.csv'
input_file = 'combined_double_links.csv'
output_file = 'converted__double_links.csv'

# read input file into a dataframe with columns names given by the first row
df = pd.read_csv(input_file, header=0)

# use the mapping in protein_to_chain to add columns next to Protein1 and Protein2 as chain1 and chain2
df['chain1'] = df['Protein1'].map(protein_to_chain)
df['chain2'] = df['Protein2'].map(protein_to_chain)
# add two columns 'Atom1' and 'Atom2' with value 'CA'
df['Atom1'] = 'CA'
df['Atom2'] = 'CA'
# drop the columns 'Protein1' and 'Protein2' and rearrange the columns in the order: chain1, Residue1, Atom1, chain2, Residue2, Atom2
df = df[['chain1', 'Residue1', 'Atom1', 'chain2', 'Residue2', 'Atom2']]
print(df.head())
# write the dataframe to the output file
df.to_csv(output_file, index=False)