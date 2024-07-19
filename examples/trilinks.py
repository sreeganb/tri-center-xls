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

import numpy as np
import os

#----------------------------------------------------------------------
# Define a bead with residue name and chain id with IMP hierarchy
# and associated degrees of freedom
#----------------------------------------------------------------------
mdl = IMP.Model()
sys = IMP.pmi.topology.System(mdl, name ='Modeling of triple crosslinking')
output = IMP.pmi.output.Output()

st1 = sys.create_state()
colors = ['red', 'blue', 'green', 'yellow', 'orange', 'purple', 'cyan', 'magenta']
seq = 'K'*1
m1 = st1.create_molecule("prot1", seq, chain_id = "B")
m1.add_representation(m1, resolutions=[1], color=colors[0])
m2 = st1.create_molecule("prot2", seq, chain_id = "C")
m2.add_representation(m2, resolutions=[1], color=colors[1])
#----------------------------------------------------------------------
# Define the path to data files
#----------------------------------------------------------------------
#pdb_dir = './data/pdb/'
pdb_dir = './data/pdb/5gjr_structure/'
fasta_dir = './data/fasta/'
xl_data = './derived_data/xl/xls_data.txt'
#----------------------------------------------------------------------
# Read in the protein structure information
# 1) Store the fasta sequences as a dictionary
#----------------------------------------------------------------------
sequences = IMP.pmi.topology.Sequences(fasta_dir + '26s_proteasome.fasta.txt')

subunit1 = st1.create_molecule("alp3", sequences["subunitalpha3"])
subunit2 = st1.create_molecule("alp6", sequences["subunitalpha6"])
subunit3 = st1.create_molecule("su10B", sequences["subunit10B"])
#subunit4 = st1.create_molecule("su6A", sequences["subunit6A"])
#subunit5 = st1.create_molecule("su8", sequences["subunit8"])

su1 = subunit1.add_structure(pdb_dir + 'alpha3.pdb', 
                             chain_id = 'n',
                             offset = 1
                            )
su2 = subunit2.add_structure(pdb_dir + 'alpha6.pdb',
                            chain_id = 'h'
                            )
su3 = subunit3.add_structure(pdb_dir + "sub10B.pdb", 
                             chain_id = 'L'
                             )
#su4 = subunit4.add_structure(pdb_dir + "subunit6A.pdb",
#                                chain_id = 'M'
#                                )
#su5 = subunit5.add_structure(pdb_dir + "subunit8.pdb",
#                                chain_id = 'J'
#                                )
#----------------------------------------------------------------------
subunit1.add_representation(su1, resolutions=[1 ,10], color=colors[2])
subunit1.add_representation(
    subunit1[:]-su1,
    # areas without structure can only be represented at one resolution
    resolutions=[1],
    # Set up spherical gaussian densities for these particles
    setup_particles_as_densities=False)
subunit2.add_representation(su2, resolutions=[1, 10], color=colors[3])
subunit2.add_representation(
    subunit2[:]-su2,
    # areas without structure can only be represented at one resolution
    resolutions=[1],
    # Set up spherical gaussian densities for these particles
    setup_particles_as_densities=False)
subunit3.add_representation(su3, resolutions=[1, 10], color=colors[4])
subunit3.add_representation(
    subunit3[:]-su3,
    # areas without structure can only be represented at one resolution
    resolutions=[1],
    # Set up spherical gaussian densities for these particles
    setup_particles_as_densities=False)
#subunit4.add_representation(su4, resolutions=[1, 10], color=colors[5])
#subunit4.add_representation(
#    subunit4[:]-su4,
#    # areas without structure can only be represented at one resolution
#    resolutions=[1],
#    # Set up spherical gaussian densities for these particles
#    setup_particles_as_densities=False)
#subunit5.add_representation(su5, resolutions=[1, 10])
#subunit5.add_representation(
#    subunit5[:]-su5,
#    # areas without structure can only be represented at one resolution
#    resolutions=[1],
#    # Set up spherical gaussian densities for these particles
#    setup_particles_as_densities=False)
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

print("non rigid parts: ", subunit1.get_non_atomic_residues())

rb1_alpha3 = dof.create_rigid_body(
    subunit1,
    #max_trans=1.0,
    #max_rot=0.5,
    nonrigid_parts=subunit1.get_non_atomic_residues())

rb2_alpha6 = dof.create_rigid_body(
    subunit2,
    #max_trans=1.0,
    #max_rot=0.5,
    nonrigid_parts=subunit2.get_non_atomic_residues())

rb3_subunit10B = dof.create_rigid_body(
    subunit3,
    #max_trans=1.0,
    #max_rot=0.5,
    nonrigid_parts=subunit3.get_non_atomic_residues())

#rb4_subunit6A = dof.create_rigid_body(
#    subunit4,
#    max_trans = 5.0,
#    max_rot = 0.5,
#    nonrigid_parts = subunit4.get_non_atomic_residues())

#rb5_subunit8 = dof.create_rigid_body(
#    subunit5,
#    max_trans = 5.0,
#    max_rot = 0.5,
#    nonrigid_parts = subunit5.get_non_atomic_residues())
#----------------------------------------------------------------------
# Connectivity restraint
#----------------------------------------------------------------------
cr = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(subunit1)
cr.add_to_model()           # add restraint to the model
output_objects.append(cr)   # add restraint to the output

cr2 = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(subunit2)
cr2.add_to_model()           # add restraint to the model
output_objects.append(cr2)   # add restraint to the output

cr3 = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(subunit3)
cr3.add_to_model()           # add restraint to the model
output_objects.append(cr3)   # add restraint to the output

#cr4 = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(subunit4)
#cr4.add_to_model()           # add restraint to the model
#output_objects.append(cr4)   # add restraint to the output

#cr5 = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(subunit5)
#cr5.add_to_model()           # add restraint to the model
#output_objects.append(cr5)   # add restraint to the output
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
                                            included_objects=[subunit1, 
                                                              subunit2, subunit3, 
                                                              m1, m2 ], #, #subunit4, subunit5, m1, m2],
                                            resolution=1000)
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
xl_weight = 4.0 

xlr = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
    root_hier=r1_hier,    # Must pass the root hierarchy to the system
    database=xldb,        # The crosslink database.
    length=11.0,          # The crosslinker plus side chain length
    resolution=1.0,       # The resolution at which to evaluate the crosslink
    slope=0.01,          # This adds a linear term to the scoring function
                          #   to bias crosslinks towards each other
    weight=xl_weight,     # Scaling factor for the restraint score.
    linker=ihm.cross_linkers.dss)  # The linker chemistry

xlr.add_to_model()
output_objects.append(xlr)
#**********************************************************************
# add some synthetic data and check if that can help our cause,
# create crosslinks between residues of the same three proteins as above
#**********************************************************************


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
                                    bounding_box=((-120,-120,-120),(120, 120, 120)),
                                    avoidcollision_rb=False)
dof.optimize_flexible_beads(100)

rex=IMP.pmi.macros.ReplicaExchange(mdl,
                                   root_hier=r1_hier,           
                                   monte_carlo_sample_objects=dof.get_movers(),
                                   replica_exchange_maximum_temperature=4.0,
                                   global_output_directory="output/",
                                   output_objects=output_objects,
                                   nframes_write_coordinates=1,
                                   monte_carlo_steps=20,
                                   number_of_frames=500,
                                   number_of_best_scoring_models=1)

rex.execute_macro()

# ihm is for depositing structures into the CIF file
# add dsso as linker for placeholder 
# Create a new database file format for crosslinks that will have an additional
# column for the dummy bead ID and this will be used by the program later to identify the dummy bead
 

#fb_alpha3 = dof.create_flexible_beads(subunit1.get_non_atomic_residues(), 
#                                      max_trans=1.0)


