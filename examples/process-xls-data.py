import pandas as pd
import numpy as np
import math
import itertools
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.inset_locator import (inset_axes, InsetPosition, mark_inset)
import math
import copy
import re
from collections import OrderedDict 
import Bio
from Bio import SeqIO
import xls_analysis_classes
from xls_analysis_classes import *
import MDAnalysis as mda
import warnings
from MDAnalysis.lib.distances import calc_bonds

#from xls_analysis_classes import *

matplotlib.rcParams.update({'font.size': 10})
d = pd.read_csv('./derived_data/xl/proteasome-tri-crosslinks.csv', sep = ',')
data = pd.DataFrame(d)
#data.insert(2, 'dummy A', 'prot1')
#data.insert(3, 'dummy XLA', 'K1')
data.insert(2, 'Subunit DB', 'prot1')
data.insert(3, 'XL DB', 'K1')
#
data.insert(6, 'dummy B', 'prot1')
data.insert(7, 'dummy XLB', 'K1')
data.insert(10, 'dummy C', 'prot1')
data.insert(11, 'dummy XLC', 'K1')
#print(data)
df1 = data.iloc[:, :4].join(data.iloc[:, 12:13]).join(data.iloc[:, 15:16]).join(data.iloc[:, 18:32]).join(data.iloc[:, 46:53])
df1.insert(6, 'Common/Gene Name B', 'prot1')
print(df1.columns)
df2 = data.iloc[:, 4:8].join(data.iloc[:, 13:14]).join(data.iloc[:, 16:17]).join(data.iloc[:, 18:25]).join(data.iloc[:, 32:39]).join(data.iloc[:, 46:53])
#print(df2.columns)
df3 = data.iloc[:, 8:12].join(data.iloc[:, 14:15]).join(data.iloc[:, 17:18]).join(data.iloc[:, 18:25]).join(data.iloc[:, 39:46]).join(data.iloc[:, 46:53])
#print(df3.columns)
df1.to_csv('proteasome-tri-crosslinks-1.csv', index = False)
df2.to_csv('proteasome-tri-crosslinks-2.csv', index = False)
df3.to_csv('proteasome-tri-crosslinks-3.csv', index = False)
print(df3)

# open a file which is down two directories and called "definitions.txt" and read it
# open this file as a pandas dataframe and name the columns as genename, genericname
heads = ['genename', 'genericname', 'chain1', 'chain2', 'subunit name']
dnames = pd.read_csv('./definitions.txt', names = heads)
#-------------------------------------------------------------------------------------
# Define a function that replaces or inserts the chain1 based on the subunit name
#-------------------------------------------------------------------------------------
def set_chain_id1(sname):
    # Check if the subunit name is NaN or empty
    if pd.isna(sname) or sname == '':
        return None
    # Find indices where the subunit name matches the sname
    indices = dnames[dnames['subunit name'].str.contains(sname, na=False)].index.tolist()
    # If there is a match, return the corresponding chain1 value
    if indices:
        return str(dnames.loc[indices[0], 'chain1'])
    return None

d['Chain ID A'] = d['Subunit A'].apply(set_chain_id1)
d['Chain ID B'] = d['Subunit B'].apply(set_chain_id1)
d['Chain ID C'] = d['Subunit C'].apply(set_chain_id1)
print(d['Chain ID B'].dropna())
columns_to_new = ['Chain ID A', 'XL A', 'Chain ID B', 'XL B']
xl_AB = d[columns_to_new].dropna()
unique_xl_AB = xl_AB.drop_duplicates()
# Split the 'XL A' and 'XL B' columns by ';' and expand into separate rows
unique_xl_AB = unique_xl_AB.assign(**{
    'XL A': unique_xl_AB['XL A'].str.split(r'[;|]'),
    'XL B': unique_xl_AB['XL B'].str.split(r'[;|]')
}).explode('XL A').explode('XL B')
unique_xl_AB['XL A'] = unique_xl_AB['XL A'].str.replace('K', '', regex=False)
unique_xl_AB['XL B'] = unique_xl_AB['XL B'].str.replace('K', '', regex=False)
# write the dataframes to csv files
unique_xl_AB.to_csv('unique_xl_AB.csv', index=False)
xl_AB.to_csv('xl_AB.csv', index=False)
columns_to_new = ['Chain ID A', 'XL A', 'Chain ID C', 'XL C']
xl_AC = d[columns_to_new].dropna()
xl_AC.to_csv('xl_AC.csv', index=False)
columns_to_new = ['Chain ID B', 'XL B', 'Chain ID C', 'XL C']
xl_BC = d[columns_to_new].dropna()
xl_BC.to_csv('xl_BC.csv', index=False)
#------------------------------------------------------------------------------------- 
# Suppress MDAnalysis warnings
warnings.filterwarnings("ignore", category=UserWarning, module="MDAnalysis")
# Load the PDB file
u = mda.Universe('./data/pdb/new-5gjr.pdb')
# Rest of the code...
# Load the file containing residue information
residue_file = './derived_data/xl/unique_xl_AB.csv'
df = pd.read_csv(residue_file)

#residues = []
#with open(residue_file, 'r') as f:
#    for line in f:
#        residue = line.strip()
#        residues.append(residue)
def calc_distance(residue_num1, residue_num2, chain_id1, chain_id2):
    # Construct selection strings for CA atoms
    selection1 = f"segid {chain_id1} and resid {residue_num1} and name CA"
    selection2 = f"segid {chain_id2} and resid {residue_num2} and name CA"

    # Select the atoms
    atom1 = u.select_atoms(selection1)
    atom2 = u.select_atoms(selection2)

    # Ensure both selections are not empty
    if len(atom1) == 0 or len(atom2) == 0:
        raise ValueError("One of the selections did not return any atoms. Please check your chain IDs and residue numbers.")

    # Get coordinates of the selected atoms
    coords1 = atom1.positions
    coords2 = atom2.positions
    # Calculate the distance between the two atoms
    print("coordinates are: ", coords1, coords2)
    print("atoms are: ", atom1, atom2)
    distance = calc_bonds(coords1, coords2)
    return distance
print(df.columns)
dis = calc_distance(df['XL A'][0], df['XL B'][0], df['Chain ID A'][0], df['Chain ID B'][0])
print("calpha-calpha distance: ", dis)

calpha_xl_distances = []
for i in range(len(df['XL A'])):
    print("what is this:" , df['XL A'][i])
    dis = calc_distance(df['XL A'][i], df['XL B'][i], df['Chain ID A'][i], df['Chain ID B'][i])
    calpha_xl_distances.append(dis)

print("calpha distances list: ", calpha_xl_distances)
    
    