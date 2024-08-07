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
m1 = st1.create_molecule("prot1", seq, chain_id = "A")
m1.add_representation(m1, resolutions=[1], color=colors[0])
m2 = st1.create_molecule("prot2", seq, chain_id = "B")
m2.add_representation(m2, resolutions=[1], color=colors[0])
m3 = st1.create_molecule("prot3", seq, chain_id = "C")
m3.add_representation(m3, resolutions=[1], color=colors[0])
m4 = st1.create_molecule("prot4", seq, chain_id = "E")
m4.add_representation(m4, resolutions=[1], color=colors[0])
m5 = st1.create_molecule("prot5", seq, chain_id = "F")
m5.add_representation(m5, resolutions=[1], color=colors[0])
#----------------------------------------------------------------------
# Define the path to data files
#----------------------------------------------------------------------
pdb_dir_1 = './data/pdb/'
pdb_dir = './data/pdb/5gjr_structure/'
fasta_dir = './data/fasta/'
xl_data = './derived_data/xl/xls_data.txt'
#----------------------------------------------------------------------
# Read in the protein structure information
# 1) Store the fasta sequences as a dictionary
# Modelling of the base subcomplex of the human proteasome
# This consists of Rpt1-6 and Rpn1, 2, 10, 13
# The structure of the base complex needs to be sequestered out of the 
# pdb or mmCIF file and used for IMP modelling
#----------------------------------------------------------------------
sequences = IMP.pmi.topology.Sequences('mod_5gjr.fasta')

#subunit1 = st1.create_molecule("alp3", sequences["subunitalpha3"])
subunit10B = st1.create_molecule("su10B", sequences["subunit10B"], chain_id = 'L')
subunit6B = st1.create_molecule("su6B", sequences["subunit6B"], chain_id = 'K')
subunit4 = st1.create_molecule("su4", sequences["subunit4"], chain_id = 'I')
subunit8 = st1.create_molecule("su8", sequences["subunit8"], chain_id = 'J')
subunit6A = st1.create_molecule("su6A", sequences["subunit6A"], chain_id = 'M')
subunit7 =st1.create_molecule("su7", sequences["subunit7"], chain_id = 'H')
subunitalpha7 = st1.create_molecule("alp7", sequences["subunitalpha7"], chain_id = 'k')
#----------------------------------------------------------------------
su7 = subunit7.add_structure(pdb_dir + "subunit7.pdb",
                            chain_id = 'H'
                            )
su10B = subunit10B.add_structure(pdb_dir + "sub10B.pdb", 
                             chain_id = 'L'
                             )
su8 = subunit8.add_structure(pdb_dir + "subunit8.pdb",
                                chain_id = 'J'
                                )
su6B = subunit6B.add_structure(pdb_dir + "subunit6B.pdb",
                                chain_id = 'K'
                                )
su4 = subunit4.add_structure(pdb_dir + "subunit4.pdb",
                                chain_id = 'I'
                                )
su6A = subunit6A.add_structure(pdb_dir + "subunit6A.pdb",
                               chain_id = 'M')
sualp7 = subunitalpha7.add_structure(pdb_dir + "alpha7.pdb",
                                    chain_id = 'k')
#----------------------------------------------------------------------
subunit10B.add_representation(su10B, resolutions=[1 ,10], color=colors[1])
subunit10B.add_representation(
    subunit10B[:]-su10B,
    # areas without structure can only be represented at one resolution
    resolutions=[1],
    # Set up spherical gaussian densities for these particles
    setup_particles_as_densities=False)
subunit8.add_representation(su8, resolutions=[1, 10], color=colors[2])
subunit8.add_representation(
    subunit8[:]-su8,
    # areas without structure can only be represented at one resolution
    resolutions=[1],
    # Set up spherical gaussian densities for these particles
    setup_particles_as_densities=False)
subunit6B.add_representation(su6B, resolutions=[1, 10], color=colors[3])
subunit6B.add_representation(
    subunit6B[:]-su6B,
    # areas without structure can only be represented at one resolution
    resolutions=[1],
    # Set up spherical gaussian densities for these particles
    setup_particles_as_densities=False)
subunit4.add_representation(su4, resolutions=[1, 10], color=colors[4])
subunit4.add_representation(
    subunit4[:]-su4,
    # areas without structure can only be represented at one resolution
    resolutions=[1],
    # Set up spherical gaussian densities for these particles
    setup_particles_as_densities=False)
subunit6A.add_representation(su6A, resolutions=[1, 10], color=colors[5])
subunit6A.add_representation(
    subunit6A[:]-su6A,
    # areas without structure can only be represented at one resolution
    resolutions=[1],
    # Set up spherical gaussian densities for these particles
    setup_particles_as_densities=False)
subunit7.add_representation(su7, resolutions=[1, 10], color=colors[6])
subunit7.add_representation(
    subunit7[:]-su7,
    # areas without structure can only be represented at one resolution
    resolutions=[1],
    # Set up spherical gaussian densities for these particles
    setup_particles_as_densities=False)
subunitalpha7.add_representation(sualp7, resolutions=[1, 10], color=colors[7])
subunitalpha7.add_representation(
    subunitalpha7[:]-sualp7,
    # areas without structure can only be represented at one resolution
    resolutions=[1],
    # Set up spherical gaussian densities for these particles
    setup_particles_as_densities=False)
#----------------------------------------------------------------------
r1_hier = sys.build()
output_objects = []
#-------------------------------------------------------------------
# Define the degrees of freedom 
#-------------------------------------------------------------------
dof = IMP.pmi.dof.DegreesOfFreedom(mdl)
#----------------------------------------------------------------------
# add the lysine bead
#----------------------------------------------------------------------
dof.create_flexible_beads(m1)
dof.create_flexible_beads(m2)
dof.create_flexible_beads(m3)
dof.create_flexible_beads(m4)
dof.create_flexible_beads(m5)

rb1_su10B = dof.create_rigid_body(
    subunit10B,
    #max_trans=1.0,
    #max_rot=0.5,
    nonrigid_parts=subunit10B.get_non_atomic_residues())

rb2_su6B = dof.create_rigid_body(
    subunit6B,
    #max_trans=1.0,
    #max_rot=0.5,
    nonrigid_parts=subunit6B.get_non_atomic_residues())

rb3_su8 = dof.create_rigid_body(
    subunit8,
    #max_trans=1.0,
    #max_rot=0.5,
    nonrigid_parts=subunit8.get_non_atomic_residues())

rb4_su4 = dof.create_rigid_body(
    subunit4,
    #max_trans = 1.0,
    #max_rot = 0.5,
    nonrigid_parts = subunit4.get_non_atomic_residues())

rb1_su6A = dof.create_rigid_body(
    subunit6A,
    nonrigid_parts=subunit6A.get_non_atomic_residues())

rb1_su7 = dof.create_rigid_body(
    subunit7,
    nonrigid_parts=subunit7.get_non_atomic_residues())
rb_sualp7 = dof.create_rigid_body(
    subunitalpha7,
    nonrigid_parts=subunitalpha7.get_non_atomic_residues())
#----------------------------------------------------------------------
# Connectivity restraint
#----------------------------------------------------------------------
def add_connectivity_restraint(subunit):
    cr = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(subunit)
    cr.add_to_model()           # add restraint to the model
    output_objects.append(cr)   # add restraint to the output
    
cr =add_connectivity_restraint(subunit10B)
add_connectivity_restraint(subunit6B)
add_connectivity_restraint(subunit8)
add_connectivity_restraint(subunit4)
add_connectivity_restraint(subunit6A)
add_connectivity_restraint(subunit7)
add_connectivity_restraint(subunitalpha7)
# -----------------------------
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
                                            included_objects=[subunit10B, subunit6B,
                                                              subunit8, subunit4, 
                                                              subunit6A, subunit7, 
                                                              subunitalpha7],
                                            resolution=1)
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
xl_weight = 120.0 

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
#----------------------------------------------------------------------
# Crosslinking based on the decomposition of the three center crosslink
# into a set of three regular crosslinks.
# Testing if this is going to provide us with better structures or not
#----------------------------------------------------------------------
#xl_data = './derived_data/xl/2_xls_data.txt'
#xldbkc = IMP.pmi.io.crosslink.CrossLinkDataBaseKeywordsConverter()
#xldbkc.set_standard_keys()
#xldb = IMP.pmi.io.crosslink.CrossLinkDataBase()
#xldb.create_set_from_file(file_name=xl_data,
#                            converter=xldbkc)
#xl_weight = 1.0

#xlr1 = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
#    root_hier=r1_hier,    # Must pass the root hierarchy to the system
#    database=xldb,          # The crosslink database.
#    length=21.0,              # The crosslinker plus side chain length
#    resolution=1,           # The resolution at which to evaluate the crosslink
#    slope=0.00001,         # This adds a linear term to the scoring function to bias crosslinks towards each other
#    weight=xl_weight,       # Scaling factor for the restraint score.
#    linker=ihm.cross_linkers.dss)  # The linker chemistry
#xlr1.add_to_model()
#output_objects.append(xlr1) 
# add a column with the ID of the dummy bead 

#####################################################
#                      SAMPLING                     #
#####################################################
# With our representation and scoring functions determined, we can now sample
# the configurations of our model with respect to the information.

# First shuffle all particles to randomize the starting point of the
# system. For larger systems, you may want to increase max_translation
sel = IMP.atom.Selection(r1_hier).get_selected_particles()
IMP.pmi.tools.shuffle_configuration(sel,
                                    max_translation=200,
                                    bounding_box=((-300,-300,-300),(300, 300, 300)),
                                    avoidcollision_rb=False)
dof.optimize_flexible_beads(100)

rex=IMP.pmi.macros.ReplicaExchange(mdl,
                                   root_hier=r1_hier,           
                                   monte_carlo_sample_objects=dof.get_movers(),
                                   replica_exchange_maximum_temperature=4.0,
                                   global_output_directory="output_new/",
                                   output_objects=output_objects,
                                   nframes_write_coordinates=1,
                                   monte_carlo_steps=25,
                                   number_of_frames=50,
                                   number_of_best_scoring_models=1)

#rex.execute_macro()

# ihm is for depositing structures into the CIF file
# add dsso as linker for placeholder 
# Create a new database file format for crosslinks that will have an additional
# column for the dummy bead ID and this will be used by the program later to identify the dummy bead
 

#fb_alpha3 = dof.create_flexible_beads(subunit1.get_non_atomic_residues(), 
#                                      max_trans=1.0)


