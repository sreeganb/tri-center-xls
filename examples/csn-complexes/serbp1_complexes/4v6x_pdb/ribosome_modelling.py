############################################
import IMP
import RMF
import IMP.atom
import IMP.rmf
import IMP.pmi
import IMP.pmi.topology
import IMP.pmi.dof
import IMP.pmi.macros
import IMP.pmi.restraints
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.restraints.crosslinking
#import IMP.pmi.restraints.pemap
from IMP.pmi.io.crosslink import CrossLinkDataBaseKeywordsConverter
#import IMP.pmi.restraints.occams

import random
import numpy as np
import glob
from sys import exit
from sys import argv
import sys
import os
#--------------------------------------------------------#
# Set up directories
#--------------------------------------------------------#
directory = os.getcwd()
data_dir = directory+'/data'
fasta_dir = data_dir + '/fasta'
pdb_dir = data_dir + '/pdb'

mdl = IMP.Model()
#--------------------------------------------------------#
# Topology
#--------------------------------------------------------#
top_file = 'topology.dat'

reader_serpb1 = IMP.pmi.topology.TopologyReader(top_file,
                                               pdb_dir = pdb_dir,
                                               fasta_dir = fasta_dir)

bs_serbp1 = IMP.pmi.macros.BuildSystem(mdl,
                                      resolutions=[1])
bs_serbp1.add_state(reader_serpb1)


hier_S1, dof_S1 = bs_serbp1.execute_macro(max_rb_trans=4.0,
                                          max_rb_rot=1.0)
mols_S1 = bs_serbp1.get_molecules()[0]

output = IMP.pmi.output.Output()
output.init_rmf("ini_all.rmf3", [hier_S1])
output.write_rmf("ini_all.rmf3")
output.close_rmf("ini_all.rmf3")

#--------------------------------------------------------#
# Represenation
#--------------------------------------------------------#
#rna1 = IMP.
#DNA0 = IMP.atom.Selection(hier_S1,
#                          resolution=1,
#                          molecule="DNA0").get_selected_particles()
##############################
# Electrostatic
##############################

#p_one = [IMP.atom.LYS, IMP.atom.ARG]
#m_one = [IMP.atom.ASP, IMP.atom.GLU, IMP.atom.DADE, IMP.atom.DCYT, IMP.atom.DGUA, IMP.atom.DTHY]
#p_five = [IMP.atom.HIS]
#
#min_distance = 8.0
#max_distance = 16.0
#
############
#ps = IMP.atom.Selection(hier_S1).get_selected_particles()
#
#for p in ps:
#    if IMP.atom.Residue(p).get_is_protein():
#        if IMP.atom.get_residue_type(p) in p_one:
#            p = IMP.atom.Charged.setup_particle(mdl, p, 1.0)
#        elif IMP.atom.get_residue_type(p) in m_one:
#            p = IMP.atom.Charged.setup_particle(mdl, p, -1.0)
#        elif IMP.atom.get_residue_type(p) in p_five:
#            p = IMP.atom.Charged.setup_particle(mdl, p, 0.5)
#        else:
#            p = IMP.atom.Charged.setup_particle(mdl, p, 0.0)
#    elif IMP.atom.Residue(p).get_is_dna():
#        p = IMP.atom.Charged.setup_particle(mdl, p, -1.0)
#        
###############################
## Lennard-Jones
###############################
#for p in ps:
#    if IMP.atom.Residue(p).get_is_protein():
#        d0 = IMP.atom.LennardJones.setup_particle(mdl, p, 0.09)
#        print(IMP.atom.Residue(p), d0.get_radius())
#    elif IMP.atom.Residue(p).get_is_dna():
#        d0 = IMP.atom.LennardJones.setup_particle(mdl, p, 0.036)
#        d0.set_radius(5.2)
#        print(IMP.atom.Residue(p), d0.get_radius())
#
################################
## Select tails
################################
#h3_tails = IMP.atom.Selection(hier_S1,
#                             resolution=1,
#                             residue_indexes=list(range(1,45)),
#                             molecule="H3").get_selected_particles()
#
#h4_tails = IMP.atom.Selection(hier_S1,
#                             resolution=1,
#                             residue_indexes=range(2,25),
#                             molecule="H4").get_selected_particles()
#
#h2a_tails = IMP.atom.Selection(hier_S1,
#                               resolution=1,
#                               residue_indexes=list(range(1,20))+list(range(120,133)),
#                               molecule="H2A").get_selected_particles()
#
#h2b_tails = IMP.atom.Selection(hier_S1,
#                               resolution=1,
#                               residue_indexes=list(range(1,40)),
#                               molecule="H2B").get_selected_particles()
#
#
#all_tails = h3_tails + h4_tails + h2a_tails + h2b_tails
#all_tails = h3_tails + h4_tails + h2a_tails + h2b_tails
#
#
################################
## Select all other
################################
#h3_core = IMP.atom.Selection(hier_S1,
#                             resolution=1,
#                             residue_indexes=list(range(45,136)),
#                             molecule="H3").get_selected_particles()
#
#h4_core = IMP.atom.Selection(hier_S1,
#                             resolution=1,
#                             residue_indexes=list(range(25,103)),
#                             molecule="H4").get_selected_particles()
#
#h2a_core = IMP.atom.Selection(hier_S1,
#                             resolution=1,
#                             residue_indexes=list(range(20,120)),
#                             molecule="H2A").get_selected_particles()
#
#h2b_core = IMP.atom.Selection(hier_S1,
#                             resolution=1,
#                             residue_indexes=list(range(40,132)),
#                             molecule="H2B").get_selected_particles()
#
#DNA0 = IMP.atom.Selection(hier_S1,
#                          resolution=1,
#                          molecule="DNA0").get_selected_particles()
#
#DNA1 = IMP.atom.Selection(hier_S1,
#                          resolution=1,
#                          molecule="DNA1").get_selected_particles()
#
#
##all_nucleosome = h3_core + h4_core + h2a_core + h2b_core + DNA0 + DNA1
#all_nucleosome =  h3_core + h4_core +  h2b_core + h2b_core + DNA0 + DNA1
#
#print('------')
#print(h2a_tails)
#print('------')
#print(h2a_core)
#
#
#def intersection(lst1, lst2): 
#    return list(set(lst1) & set(lst2))
#
#print(intersection(h3_tails, h3_core))
#print(intersection(h4_tails, h4_core))
#print(intersection(h2a_tails, h2a_core))
#print(intersection(h2b_tails, h2b_core))

###############################
# Connectivity
###############################

#output_objects = [] # keep a list of functions that need to be reported
#sample_objects = []
#rmf_restraints = []
#display_restraints = []
#
#
#crs = []
#for molname in mols_S1:
#    for mol in mols_S1[molname]:
#        copy_n = IMP.atom.Copy(mol.get_hierarchy()).get_copy_index()
#        cr = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(mol)
#        cr.set_label(mol.get_name()+'.'+str(copy_n))
#        cr.add_to_model()
#        output_objects.append(cr)
#        crs.append(cr)
#        
#    
###############################
## Shuffle
###############################    
#IMP.pmi.tools.shuffle_configuration(hier_S1)
#                                    #max_translation=200,
#                                    #bounding_box=((-80, -80, -80), (80, 80, 80)))
#dof_S1.optimize_flexible_beads(200)
#
#
#print(dof_S1.get_movers())
#
############################## SAMPLING ##############################
## Run replica exchange Monte Carlo sampling
##taskid = int(sys.argv[1])
#mc1 = IMP.pmi.macros.ReplicaExchange0(mdl,
#                                    root_hier=hier_S1,                           
#                                    crosslink_restraints=rmf_restraints,          
#                                    monte_carlo_sample_objects=dof_S1.get_movers()[1:],   
#                                    replica_exchange_maximum_temperature=5.0,
#                                    global_output_directory='output',
#                                    output_objects=output_objects,
#                                    monte_carlo_steps=20,
#                                    number_of_frames=50,
#                                    number_of_best_scoring_models=0)
#
#mc1.execute_macro()
#rex1 = mc1.get_replica_exchange_object()
#
###############################
## Angle restraint
###############################
#rar = IMP.pmi.restraints.stereochemistry.ResidueAngleRestraint(
#    objects=all_tails,
#    strength=1.0)
#rar.set_weight(w_A)
#rar.add_to_model()
#output_objects.append(rar)
#
###############################
## Debye-Huckel
###############################
#
#dh = DebyeHuckel(
#    included_objects=all_tails,
#    resolution=1,
#    debye_length=7.84,
#    dielectric_constant=80.0
#)
#
#dh.add_to_model()
#dh.set_label('DH_intra')
#dh.set_weight(w_DH)
#output_objects.append(dh)
#print(dh.evaluate())
#
## bipartite
#dhb = DebyeHuckel(
#    included_objects = all_tails,
#    other_objects = all_nucleosome,
#    resolution=1,
#    debye_length=7.84,
#    dielectric_constant=80.0
#)
#
#dhb.add_to_model()
#dhb.set_label('DH_inter')
#dhb.set_weight(1.0*w_DH)
#output_objects.append(dhb)
#print(dhb.evaluate())
#
###############################
## Lennard-Jones
###############################
#
#LJ = LennardJonesExcludedVolume(
#    included_objects = all_tails,
#    resolution=1
#)
#
#LJ.add_to_model()
#LJ.set_label('LJ_intra')
#LJ.set_weight(w_LJ)
#output_objects.append(LJ)
#
#LJb = LennardJonesExcludedVolume(
#    included_objects = all_tails,
#    other_objects = all_nucleosome,
#    resolution=1
#)
#
#LJb.add_to_model()
#LJb.set_label('LJ_inter')
#LJb.set_weight(w_LJ)
#output_objects.append(LJb)
#
###############################
## pEMAP
###############################
#if include_pemap:
#    pemap = IMP.pmi.restraints.pemap.pEMapRestraint(hier_S1,
#                                                    mic_file,
#                                                    sigma_init=5.0)
#
#    
#    pemap.add_to_model()
#    pemap.set_weight(1.0)
#    output_objects.append(pemap)
#    print('output pemap: ', pemap.get_output())
#    dof_S1.get_nuisances_from_restraint(pemap)
#
###############################
## Sampling
###############################
#
#mc2 = IMP.pmi.macros.ReplicaExchange0(mdl,
#                                      root_hier=hier_S1,
#                                      crosslink_restraints=rmf_restraints,
#                                      monte_carlo_sample_objects=dof_S1.get_movers()[1:],
#                                      replica_exchange_maximum_temperature=3.0,
#                                      global_output_directory="output_A%s_LG%s_DH%s/"%(w_A,w_LJ,w_DH),
#                                      output_objects=output_objects,
#                                      monte_carlo_steps=10,
#                                      number_of_frames=100000,
#                                      number_of_best_scoring_models=0,
#                                      replica_exchange_object = rex1)
#
#mc2.execute_macro()

exit()


















