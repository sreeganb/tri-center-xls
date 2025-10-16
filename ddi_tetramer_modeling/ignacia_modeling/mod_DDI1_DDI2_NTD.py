############################################
# Modeling ddi2/ddi2 dimer and tetramer
# occams restraint symmetry
#
# iecheverria - Sali Lab - UCSF
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
import IMP.pmi.restraints.basic
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.restraints.crosslinking
#import IMP.pmi.restraints.occams
from IMP.pmi.io.crosslink import CrossLinkDataBaseKeywordsConverter

import random
import numpy as np
import glob
import os
from sys import exit
from sys import argv

#from ambiguous_xls_restraint_new import *
###################### SYSTEM SETUP #####################

#top_dir = '/wynton/home/sali/ignacia/ddi_ambiguity/modeling/mod_ddi2_ddi2_tetramer'
top_dir = './'

mdl = IMP.Model()

#cldbkc=CrossLinkDataBaseKeywordsConverter()
#cldbkc.set_protein1_key("Prot A")
#cldbkc.set_protein2_key("Prot B")
#cldbkc.set_residue1_key("Residue A")
#cldbkc.set_residue2_key("Residue B")
#cldbkc.set_unique_id_key("Iid")
#cldbkc.set_psi_key("Score")

cldbkc=CrossLinkDataBaseKeywordsConverter()
cldbkc.set_protein1_key("Protein1")
cldbkc.set_protein2_key("Protein2")
cldbkc.set_residue1_key("Residue1")
cldbkc.set_residue2_key("Residue2")
cldbkc.set_unique_id_key("Iid")
cldbkc.set_psi_key("Score")

cldb=IMP.pmi.io.crosslink.CrossLinkDataBase(cldbkc)

def create_hier_sys(top_file, name = 'sys'):
    topology = os.path.join(top_dir,top_file)
    reader = IMP.pmi.topology.TopologyReader(topology,
                                             pdb_dir = os.path.join(top_dir,'data'),
                                             fasta_dir = os.path.join(top_dir,'data'))

    
    bs = IMP.pmi.macros.BuildSystem(mdl,
                                    resolutions=[1])
    bs.add_state(reader,
		 keep_chain_id=True)

    hier, dof = bs.execute_macro()
    mols = bs.get_molecules()[0]

    output = IMP.pmi.output.Output()
    output.init_rmf(f"ini_{name}.rmf3", [hier])
    output.write_rmf(f"ini_{name}.rmf3")
    output.close_rmf(f"ini_{name}.rmf3")

    return hier, dof, mols

def add_connectivity(mols, name = '1'):
    for molname in mols:
        for mol in mols[molname]:
            copy_n = IMP.atom.Copy(mol.get_hierarchy()).get_copy_index()
            cr = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(mol)
            cr.set_label(name+'.'+mol.get_name()+'.'+str(copy_n))
            cr.add_to_model()
            output_objects.append(cr)
            crs.append(cr)

def add_excluded_volume(mols,weight = 1., name = '1'):
    evr = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(included_objects=mols.values(),
                                                                  resolution=10)
    evr.set_label(f'{name}')
    evr.add_to_model()
    evr.set_weight(weight)
    output_objects.append(evr)

def add_crosslinks(file_in,hier, dof, weight = 1., length = 21.0, name = '1'):
    cldb=IMP.pmi.io.crosslink.CrossLinkDataBase(cldbkc)
    cldb.create_set_from_file(file_in)

    xls = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(root_hier=hier,
                                                                                database=cldb,
                                                                                resolution=1,
                                                                                length=length,
                                                                                label=name)

    xls.add_to_model()
    xls.set_weight(weight)
    rmf_restraints.append(xls)
    output_objects.append(xls)
    dof.get_nuisances_from_restraint(xls)

def add_amb_crosslinks(file_in, hiers, dof, length=21.0, weight = 1., slope = 0.001, name = '1'):


    
    xls = XLs_amb_Restraint(hiers,
                            file_in,
                            sigma_init = 3.0,
			    slope = slope,
                            length = length,
                            label= name)
                                                                                

    xls.add_to_model()
    xls.set_weight(weight)
    rmf_restraints.append(xls)
    output_objects.append(xls)
    dof.get_nuisances_from_restraint(xls)

    return xls

def add_barrier(hier, name = '1'):
  
    eb = IMP.pmi.restraints.basic.ExternalBarrier(hierarchies=hier, radius=200)
    #eb.add_label(f'{name}')
    eb.add_to_model()


def add_Symmetry(hier, dof, selection, weight=1., name = '1'):
  
    sy = SymmetryRestraint(hier=hier,
                           resolution=5,
                           selection=selection,
                           label=name)

    sy.add_to_model()
    sy.set_weight(weight)
    rmf_restraints.append(sy)
    output_objects.append(sy)
    dof.get_nuisances_from_restraint(sy)

def add_COM(hier, dof, selection1, selection2,  distance=30., weight=1., name = '1'):
  
    com = COMRestraint(hier=hier,
                       selection1=selection1,
                       selection2=selection2,
                       distance=distance,
                       label=name)

    com.add_to_model()
    com.set_weight(weight)
    rmf_restraints.append(com)
    output_objects.append(com)
    dof.get_nuisances_from_restraint(com)
    
###############################
# All species 
###############################
hier_S1, dof_S1, mols_S1 = create_hier_sys(f'{top_dir}/top_DDI1_NTD.dat', name = 'tetramer1')
hier_S2, dof_S2, mols_S2 = create_hier_sys(f'{top_dir}/top_DDI2_NTD.dat', name = 'tetramer2')

##############################
# Combined hierarchy
##############################

p = IMP.Particle(mdl)
hier_all = IMP.atom.Hierarchy.setup_particle(p)
hier_all.add_child(hier_S1)
hier_all.add_child(hier_S2)
hier_all.set_name('System')

states = IMP.atom.get_by_type(hier_all,IMP.atom.STATE_TYPE)
print('All states:', states)

##############################
# Connectivity
##############################
output_objects = [] 
sample_objects = []
rmf_restraints = []

crs = []

add_connectivity(mols_S1, name='tetramer1')
add_connectivity(mols_S2, name='tetramer2')

##############################
# Excluded Volume
##############################
add_excluded_volume(mols_S1, name='tetramer1')
add_excluded_volume(mols_S2, name='tetramer2')

##############################
# External barrier
##############################

add_barrier(hier_S1, name = 'tetramer1')
add_barrier(hier_S2, name = 'tetramer2')

##############################
# XLs
##############################

#xls_tetramer = f'{top_dir}/data/mixed_DDI1_dimer.csv'
#add_amb_crosslinks(xls_tetramer, [hier_S1], dof_S1, weight = 5.0, length = 21., slope=0.01 ,name='tetramer1')
#xls_tetramer = f'{top_dir}/data/mixed_DDI2_dimer.csv'
#add_amb_crosslinks(xls_tetramer, [hier_S2], dof_S2, weight = 5.0, length = 21., slope=0.01 ,name='tetramer2')

##############################
# Write everything
##############################  
    
output = IMP.pmi.output.Output()
output.init_rmf("ini_all.rmf3", [hier_all])
output.write_rmf("ini_all.rmf3")
output.close_rmf("ini_all.rmf3")

##############################
# Shuffle
##############################    
IMP.pmi.tools.shuffle_configuration(hier_S1, max_translation=100)
IMP.pmi.tools.shuffle_configuration(hier_S2, max_translation=100)

movers_all = dof_S1.get_movers() + dof_S2.get_movers()

##############################
# MC
##############################    

rex=IMP.pmi.macros.ReplicaExchange(mdl,
                                   root_hier=hier_all,                                        
                                   monte_carlo_sample_objects=movers_all,   
                                   replica_exchange_maximum_temperature=4.0,
                                   global_output_directory="output/",
                                   output_objects=output_objects,
                                   monte_carlo_steps=20,
                                   number_of_frames=10000,
                                   number_of_best_scoring_models=0)

rex.execute_macro()

exit()

