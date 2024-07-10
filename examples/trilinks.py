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

seq = 'K'*1
m1 = st1.create_molecule("prot1", seq, chain_id = "B")
m1.add_representation(m1, resolutions=[1])
#----------------------------------------------------------------------
# Define the path to data files
#----------------------------------------------------------------------
pdb_dir = './data/pdb/'
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

su1 = subunit1.add_structure(pdb_dir + 'alpha3.pdb', 
                             chain_id = 'M',
                             #res_range = (81, 433),
                             offset = 1
                            )
su2 = subunit2.add_structure(pdb_dir + 'alpha6.pdb',
                            chain_id = 'G'
                            #res_range = (15, 389),
                            #offset = 92
                            )
su3 = subunit3.add_structure(pdb_dir + "26s_protease_subunit10B.pdb", 
                             chain_id = 'E'
                             #res_range = (39, 417), 
                             #offset = -38
                             )

#----------------------------------------------------------------------
subunit1.add_representation(su1, resolutions=[1 ,10])
subunit1.add_representation(
    subunit1[:]-su1,
    # areas without structure can only be represented at one resolution
    resolutions=[1],
    # Set up spherical gaussian densities for these particles
    setup_particles_as_densities=False)
subunit2.add_representation(su2, resolutions=[1, 10])
subunit2.add_representation(
    subunit2[:]-su2,
    # areas without structure can only be represented at one resolution
    resolutions=[1],
    # Set up spherical gaussian densities for these particles
    setup_particles_as_densities=False)
subunit3.add_representation(su3, resolutions=[1, 10])
subunit3.add_representation(
    subunit3[:]-su3,
    # areas without structure can only be represented at one resolution
    resolutions=[1],
    # Set up spherical gaussian densities for these particles
    setup_particles_as_densities=False)

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

print("non rigid parts: ", subunit1.get_non_atomic_residues())

rb1_alpha3 = dof.create_rigid_body(
    subunit1,
    max_trans=1.0,
    max_rot=0.5,
    nonrigid_parts=subunit1.get_non_atomic_residues())

rb2_alpha6 = dof.create_rigid_body(
    subunit2,
    max_trans=1.0,
    max_rot=0.5,
    nonrigid_parts=subunit2.get_non_atomic_residues())
rb3_subunit10B = dof.create_rigid_body(
    subunit3,
    max_trans=1.0,
    max_rot=0.5,
    nonrigid_parts=subunit3.get_non_atomic_residues())

cr = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(subunit1)
cr.add_to_model()           # add restraint to the model
output_objects.append(cr)   # add restraint to the output

cr2 = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(subunit2)
cr2.add_to_model()           # add restraint to the model
output_objects.append(cr2)   # add restraint to the output

cr3 = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(subunit3)
cr3.add_to_model()           # add restraint to the model
output_objects.append(cr3)   # add restraint to the output

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
                                            included_objects=[subunit1, subunit2, subunit3],
                                            resolution=1000)
output_objects.append(evr)

xldbkc = IMP.pmi.io.crosslink.CrossLinkDataBaseKeywordsConverter()
xldbkc.set_standard_keys()
xldb = IMP.pmi.io.crosslink.CrossLinkDataBase()
xldb.create_set_from_file(file_name=xl_data,
                          converter=xldbkc)
xl_weight = 1.0 

xlr = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
    root_hier=r1_hier,    # Must pass the root hierarchy to the system
    database=xldb,          # The crosslink database.
    length=25,              # The crosslinker plus side chain length
    resolution=1,           # The resolution at which to evaluate the crosslink
    slope=0.000001,         # This adds a linear term to the scoring function
                            #   to bias crosslinks towards each other
    weight=xl_weight,       # Scaling factor for the restraint score.
    linker=ihm.cross_linkers.dss)  # The linker chemistry

xlr.add_to_model()
output_objects.append(xlr)

# add a column with the ID of the dummy bead 

#####################################################
#                      SAMPLING                     #
#####################################################
# With our representation and scoring functions determined, we can now sample
# the configurations of our model with respect to the information.

# First shuffle all particles to randomize the starting point of the
# system. For larger systems, you may want to increase max_translation
IMP.pmi.tools.shuffle_configuration(r1_hier,
                                    max_translation=50)

# Shuffling randomizes the bead positions. It's good to
# allow these to optimize first to relax large connectivity
# restraint scores.  100-500 steps is generally sufficient.
dof.optimize_flexible_beads(500)

evr.add_to_model()

# Run replica exchange Monte Carlo sampling
rex = IMP.pmi.macros.ReplicaExchange(
    mdl,
    # pass the root hierarchy
    root_hier=r1_hier,
    # pass all objects to be moved ( almost always dof.get_movers() )
    monte_carlo_sample_objects=dof.get_movers(),
    # The output directory for this sampling run.
    global_output_directory='run_manual1/output/',
    # Items in output_objects write information to the stat file.
    output_objects=output_objects,
    # Number of MC steps between writing frames
    monte_carlo_steps=10,
    # set >0 to store best PDB files (but this is slow)
    number_of_best_scoring_models=0,
    # Total number of frames to run / write to the RMF file.
    number_of_frames=10000)

# Ok, now we finally do the sampling!
rex.execute_macro()


# ihm is for depositing structures into the CIF file
# add dsso as linker for placeholder 

#fb_alpha3 = dof.create_flexible_beads(subunit1.get_non_atomic_residues(), 
#                                      max_trans=1.0)


