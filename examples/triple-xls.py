#!/usr/bin/env python

import IMP
import IMP.isd
import IMP.core
import IMP.pmi.restraints.basic
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.restraints.crosslinking
import IMP.pmi.tools
import IMP.pmi.topology
import IMP.pmi.dof
import ihm.cross_linkers
import IMP.pmi.macros

from Bio import PDB, SeqIO
import warnings

import numpy as np
import pandas as pd
import os
import re

#----------------------------------------------------------------------
# Define a bead with residue name and chain id with IMP hierarchy
# and associated degrees of freedom
#----------------------------------------------------------------------
mdl = IMP.Model()
sys = IMP.pmi.topology.System(mdl, name ='Modeling of triple crosslinking')
output = IMP.pmi.output.Output()

st1 = sys.create_state()
colors = ['red', 'blue', 'green', 'yellow', 'orange', 
          'purple', 'cyan', 'magenta', 'black', 
            'pink', 'brown', 'gold']
seq = 'K'*1
df_xl = pd.read_csv('./derived_data/xl/final_triple_xls.csv')

nbeads = int(len(df_xl['Protein2'])/3) 
print("number of beads: ", nbeads)

for i in range(nbeads):
    m = st1.create_molecule(f"prot{i+1}", seq, chain_id = chr(ord('A')+i))
    m.add_representation(m, resolutions=[1], color=colors[i])