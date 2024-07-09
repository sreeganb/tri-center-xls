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

import numpy as np
import os

#----------------------------------------------------------------------
# Define a bead with residue name and chain id with IMP hierarchy
# and associated degrees of freedom
#----------------------------------------------------------------------
mdl = IMP.Model()
seq = 'K'*1
output = IMP.pmi.output.Output()
s1 = IMP.pmi.topology.System(mdl)
st1 = s1.create_state()
m1 = st1.create_molecule("prot1", seq, chain_id = "B")
m1.add_representation(m1, resolutions=[1])
r1_hier = s1.build()
sys = IMP.pmi.topology.System(mdl, name ='Modeling of triple crosslinking')
dof_s1 = IMP.pmi.dof.DegreesOfFreedom(mdl)
dof_s1.create_flexible_beads(m1)
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
#sequences = IMP.pmi.topology.Sequences(fasta_dir + '26s_proteasome_6B_7_10B.fasta.txt')
seq1 = IMP.pmi.topology.Sequences(fasta_dir + '26s_protease_subunit6B.fasta.txt')
seq2 = IMP.pmi.topology.Sequences(fasta_dir + '26s_protease_subunit7.fasta.txt')
seq3 = IMP.pmi.topology.Sequences(fasta_dir + '26s_protease_subunit10B.fasta.txt')

subunit1 = st1.create_molecule("D", seq1["subunit6B"])
subunit2 = st1.create_molecule("A", seq2["subunit7"])
subunit3 = st1.create_molecule("E", seq3["subunit10B"])

su1 = subunit1.add_structure(pdb_dir + "26s_protease_subunit6B.pdb", 
                             chain_id = "D",
                             #res_range = (39, 417), 
                             #offset = -38
                             )
su2 = subunit2.add_structure(pdb_dir + '26s_protease_subunit7.pdb', 
                             chain_id = 'A',
                             #res_range = (81, 433),
                             #offset = -80
                            )
su3 = subunit3.add_structure(pdb_dir + '26s_protease_subunit10B.pdb',
                            chain_id = 'E',
                            #res_range = (15, 389),
                            #offset = -14
                            )
#----------------------------------------------------------------------
subunit1.add_representation(su1, resolutions=[1 ,10],
                            density_residues_per_component=10
                            )
subunit1.add_representation(
    subunit1[:]-su1,
    # areas without structure can only be represented at one resolution
    resolutions=[1],
    # Set up spherical gaussian densities for these particles
    setup_particles_as_densities=True)
subunit2.add_representation(su2, resolutions=[1, 10],
                            density_residues_per_component=10)
subunit2.add_representation(
    subunit2[:]-su2,
    # areas without structure can only be represented at one resolution
    resolutions=[1],
    # Set up spherical gaussian densities for these particles
    setup_particles_as_densities=True)
subunit3.add_representation(su3, resolutions=[1, 10],
                            density_residues_per_component=10)
subunit3.add_representation(
    subunit3[:]-su3,
    # areas without structure can only be represented at one resolution
    resolutions=[1],
    # Set up spherical gaussian densities for these particles
    setup_particles_as_densities=True)

r1_hier = sys.build()