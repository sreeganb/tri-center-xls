#!/usr/bin/env python3.11

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
          'purple', 'cyan', 'magenta', 
          'pink', 'brown', 'gold']
seq = 'K'*1
df_xl = pd.read_csv('./derived_data/xl/final_triple_xls.csv')

nbeads = int(len(df_xl['Protein2'])/3) # number of lysine residues
nbeads = 13 # hardcoded for now
#----------------------------------------------------------------------
# Create a bead for each of the triple crosslinks in the experiment
#----------------------------------------------------------------------
def create_molecule(state, name, seq, chain_id, color_index):
    """Create a molecule and add it to the system with the specified color."""
    molecule = state.create_molecule(name, seq, chain_id=chain_id)
    molecule.add_representation(molecule, resolutions=[1], color=colors[color_index])
    return molecule
# Create molecules 
molecules = []
for i in range(nbeads):
  molecule = create_molecule(st1, f"p{i+1}", 'K'*1, chr(65+i), 10)
  molecules.append(molecule)
#----------------------------------------------------------------------
# Define the path to data files
#----------------------------------------------------------------------
pdb_dir = './data/pdb/'
fasta_dir = './data/fasta/'
xl_data = './derived_data/xl/final_triple_xls.csv'
#----------------------------------------------------------------------
# Read in the protein structure information
# create a list of protein names that are part of the base of the 
# proteasome and that will be modelled here, I have already extracted the 
# right chain IDs so I will use them here (hardcoded but we need to get it
# working first, so its fine)
#----------------------------------------------------------------------
sequences = IMP.pmi.topology.Sequences('mod_5gjr.fasta')
def create_and_process_molecules(proteins, chain_ids, sequences, st, colors):
  """
  Creates and processes molecules based on given protein names, chain IDs, and other parameters.

  Args:
    proteins: A list of protein names.
    chain_ids: A list of chain IDs corresponding to the proteins.
    sequences: An IMP.pmi.topology.Sequences object.
    st: An IMP state.
    pdb_dir: The directory containing PDB files.
    colors: A list of colors.

  Returns:
    A list of created molecules.
  """

  molecules = []
  sub_names = []
  pdb_file = "./data/pdb/base_proteasome.pdb"
  for protein, chain_id, color_idx in zip(proteins, chain_ids, range(len(colors))):
    #print(f"Creating molecule for {protein} with chain ID {chain_id} and color {colors[color_idx]}")
    molecule_name = f"{protein}"  # Assuming protein names end with a number
    print(f"Creating molecule {molecule_name}")
    molecule = st.create_molecule(molecule_name, sequences[protein], chain_id=chain_id)
    molecules.append(molecule)
    sub_names.append(molecule_name)

    #pdb_file = f"{pdb_dir}{protein}.pdb"
    structure = molecule.add_structure(pdb_file, chain_id=chain_id)

    molecule.add_representation(structure, resolutions=[1, 10], color=colors[color_idx])
    molecule.add_representation(
        molecule[:] - structure,
        resolutions=[50],
        setup_particles_as_densities=False
    )

  return molecules, sub_names

proteins = ['Rpt1', 'Rpt2', 'Rpt3', 'Rpt4', 'Rpt5', 'Rpt6', 'Rpn2']
chain_ids = ['v', 'w', 'y', 'z', '0', 'x', '1']

protein_subunits, sub_names = create_and_process_molecules(proteins, chain_ids, sequences, st1, colors)
molecules.extend(protein_subunits)
print(molecules)
#----------------------------------------------------------------------
r1_hier = sys.build()
output_objects = []
#-------------------------------------------------------------------
# Define the degrees of freedom 
# create flexible beads for all the nbeads number of lysine residues
# create a rigid body for each of the protein subunits with names rb_1, rb_2, ...
#-------------------------------------------------------------------
dof = IMP.pmi.dof.DegreesOfFreedom(mdl)
# Create flexible beads for all molecules
for molecule in molecules[:nbeads]:
    dof.create_flexible_beads(molecule)

# Create rigid bodies for all subunits
for subunit, mol in zip(sub_names, molecules[nbeads:]):
    rb_name = f"rb_{subunit}"  # Create a dynamic name for the rigid body
    locals()[rb_name] = dof.create_rigid_body(
        mol,
        nonrigid_parts=mol.get_non_atomic_residues()
    )
#----------------------------------------------------------------------
# Add the connectivity restraints
#----------------------------------------------------------------------
def add_connectivity_restraints(subunits):
  """Adds connectivity restraints for a list of subunits.

  Args:
    subunits: A list of subunits.
  """

  for subunit in subunits:
    cr = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(subunit)
    cr.add_to_model()
    output_objects.append(cr)

# Assuming you have a list of subunits:
subunits = molecules[nbeads:]
add_connectivity_restraints(subunits)
# %%%%% EXCLUDED VOLUME RESTRAINT
#
# Keeps particles from occupying the same area in space.
# Here, we pass a list of both molecule chains to included_objects to
# apply this to every residue.
# We could also have passed root_hier to obtain the same behavior.
#
# resolution=1000 applies this expensive restraint to the lowest
# resolution for each particle.
evr = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(
                                            included_objects=molecules[nbeads:],
                                            resolution=5)
evr.add_to_model()
output_objects.append(evr)
#----------------------------------------------------------------------
# Crosslinking restraint
#----------------------------------------------------------------------
xldbkc = IMP.pmi.io.crosslink.CrossLinkDataBaseKeywordsConverter()
xldbkc.set_standard_keys()
xldb = IMP.pmi.io.crosslink.CrossLinkDataBase()
xldb.create_set_from_file(file_name=xl_data,
                          converter=xldbkc)
xl_weight = 130.0 
xlr = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
    root_hier=r1_hier,    # Must pass the root hierarchy to the system
    database=xldb,        # The crosslink database.
    length=21.0,          # The crosslinker plus side chain length
    resolution=1.0,       # The resolution at which to evaluate the crosslink
    slope=0.01,          # This adds a linear term to the scoring function
                          #   to bias crosslinks towards each other
    weight=xl_weight,     # Scaling factor for the restraint score.
    linker=ihm.cross_linkers.dss)  # The linker chemistry

xlr.add_to_model()
output_objects.append(xlr)
#####################################################
#                      SAMPLING                     #
#####################################################
# With our representation and scoring functions determined, we can now sample
# the configurations of our model with respect to the information.

# First shuffle all particles to randomize the starting point of the
# system. For larger systems, you may want to increase max_translation
sel = IMP.atom.Selection(r1_hier).get_selected_particles()
IMP.pmi.tools.shuffle_configuration(sel,
                                    max_translation=300,
                                    bounding_box=((-400,-400,-400),(400, 400, 400)),
                                    avoidcollision_rb=False)
dof.optimize_flexible_beads(100)
rex=IMP.pmi.macros.ReplicaExchange(mdl,
                                   root_hier=r1_hier,           
                                   monte_carlo_sample_objects=dof.get_movers(),
                                   replica_exchange_maximum_temperature=4.0,
                                   global_output_directory="output_new/",
                                   output_objects=output_objects,
                                   nframes_write_coordinates=1,
                                   monte_carlo_steps=20,
                                   number_of_frames=1000,
                                   number_of_best_scoring_models=1)

rex.execute_macro()