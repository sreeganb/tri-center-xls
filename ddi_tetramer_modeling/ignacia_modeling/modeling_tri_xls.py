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
import pandas as pd

from ambiguous_xls_restraint_new import *
###################### SYSTEM SETUP #####################

#top_dir = '/wynton/home/sali/ignacia/ddi_ambiguity/modeling/mod_ddi2_ddi2_tetramer'
top_dir = './'

mdl = IMP.Model()

#=============================
# trifunctional crosslinker database
#============================
dir = os.getcwd()
xl_dir = os.path.join(dir,'data')
xl_fil = os.path.join(xl_dir,'reduced_ddi_trifunctional.csv')
#=============================

#cldbkc=CrossLinkDataBaseKeywordsConverter()
#cldbkc.set_protein1_key("Prot A")
#cldbkc.set_protein2_key("Prot B")
#cldbkc.set_residue1_key("Residue A")
#cldbkc.set_residue2_key("Residue B")
#cldbkc.set_unique_id_key("Iid")
#cldbkc.set_psi_key("Score")
#
#cldb=IMP.pmi.io.crosslink.CrossLinkDataBase(cldbkc)

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

#----------------------------------------------------------------------
# Trifunctional crosslink restraints with ambiguous copy assignments
#----------------------------------------------------------------------
def add_distance_restraints_from_csv(file_in, hier, min_distance=0.0, max_distance=19.0, 
                                     resolution=1.0, kappa=1.0, weight=1.0, name='1'):
    """
    Add trifunctional XL restraints handling copy number ambiguity.
    
    Logic (matching boss's bifunctional XL approach):
    - Same protein, same residue = MUST be different copies (no self-crosslink)
    - Same protein, different residues = Could be intra OR inter-copy (ambiguous)
    - Different proteins = All copy combinations possible
    
    For each trifunctional XL, creates 3 distance restraints (forming a triangle)
    for ALL valid copy number combinations.
    """
    import itertools
    
    xldata = pd.read_csv(file_in, header=0)
    total_restraints = 0
    
    for idx, row in xldata.iterrows():
        # Extract protein names and residue numbers
        proteins = [row['Protein1'], row['Protein2'], row['Protein3']]
        residues = [int(''.join(filter(str.isdigit, str(row[f'Residue{i}'])))) 
                   for i in [1, 2, 3]]
        
        # Determine valid copy number combinations for this trifunctional XL
        valid_combos = _get_valid_copy_combinations(proteins, residues)
        
        print(f"\nXL {idx}: {proteins[0]}:{residues[0]} - {proteins[1]}:{residues[1]} - {proteins[2]}:{residues[2]}")
        print(f"  Valid copy combinations: {len(valid_combos)}")
        
        # Add restraints for each valid combination
        for (copy1, copy2, copy3) in valid_combos:
            # Get particles for this specific copy combination
            particles = _select_particles(hier, proteins, residues, 
                                         [copy1, copy2, copy3], resolution)
            
            if particles is None:  # Skip if any particle not found
                continue
            
            # Create 3 pairwise distance restraints (triangle constraint)
            pairs = [(0,1), (1,2), (0,2)]  # Indices into particles list
            pair_names = ['1-2', '2-3', '1-3']
            
            for (i, j), pair_name in zip(pairs, pair_names):
                tup_i = [residues[i], residues[i], proteins[i], [copy1, copy2, copy3][i]]
                tup_j = [residues[j], residues[j], proteins[j], [copy1, copy2, copy3][j]]
                
                dr = IMP.pmi.restraints.basic.DistanceRestraint(
                    hier, tup_i, tup_j, min_distance, max_distance, 
                    resolution=resolution, kappa=kappa)
                
                label = f"{name}.{proteins[i]}.{[copy1,copy2,copy3][i]}:{residues[i]}-" \
                       f"{proteins[j]}.{[copy1,copy2,copy3][j]}:{residues[j]}"
                dr.set_label(label)
                dr.add_to_model()
                dr.set_weight(weight)
                output_objects.append(dr)
                total_restraints += 1
            
            print(f"  Added combo: {proteins[0]}.{copy1}:{residues[0]} - "
                  f"{proteins[1]}.{copy2}:{residues[1]} - {proteins[2]}.{copy3}:{residues[2]}")
    
    print(f"\n=== Total: {total_restraints} distance restraints from {len(xldata)} trifunctional XLs ===\n")


def _get_valid_copy_combinations(proteins, residues):
    """
    Determine all valid (copy1, copy2, copy3) combinations.
    
    Rules:
    1. Same protein + same residue = different copies required
    2. Same protein + different residues = any copy combination (ambiguous)
    3. Different proteins = any copy combination
    """
    def copy_options(prot1, res1, prot2, res2):
        """Get possible copy pairs for a protein-residue pair"""
        if prot1 != prot2:
            return [(0,0), (0,1), (1,0), (1,1)]  # Different proteins: all combos
        elif res1 == res2:
            return [(0,1), (1,0)]  # Same residue: MUST be different copies
        else:
            return [(0,0), (1,1), (0,1), (1,0)]  # Different residues: ambiguous
    
    # Get possible copy pairs for each edge of the triangle
    combos_12 = copy_options(proteins[0], residues[0], proteins[1], residues[1])
    combos_23 = copy_options(proteins[1], residues[1], proteins[2], residues[2])
    combos_13 = copy_options(proteins[0], residues[0], proteins[2], residues[2])
    
    # Find consistent (copy1, copy2, copy3) combinations
    # Consistency: copy numbers must match when same position appears in multiple pairs
    valid = set()
    for (c1_12, c2_12) in combos_12:
        for (c2_23, c3_23) in combos_23:
            for (c1_13, c3_13) in combos_13:
                # Check if copy assignments are consistent across all three pairs
                if c1_12 == c1_13 and c2_12 == c2_23 and c3_23 == c3_13:
                    valid.add((c1_12, c2_12, c3_23))
    
    return sorted(valid)


def _select_particles(hier, proteins, residues, copies, resolution):
    """
    Select particles for given proteins/residues/copies.
    Returns None if any particle missing (with warning).
    """
    particles = []
    for prot, res, copy in zip(proteins, residues, copies):
        sel = IMP.atom.Selection(hier, state_index=0, molecule=prot, 
                                residue_index=res, copy_index=copy, 
                                resolution=resolution)
        p = sel.get_selected_particles()
        
        if not p:
            print(f"  Warning: No particles for {prot}.{copy}:{res} - skipping this combination")
            return None
        
        particles.append(p[0])
    
    return particles
#----------------------------------------------------------------------

#----------------------------------------------------------------------
###############################
# All species 
###############################
#hier_S1, dof_S1, mols_S1 = create_hier_sys(f'{top_dir}/top_DDI1_NTD.dat', name = 'tetramer1')
#hier_S2, dof_S2, mols_S2 = create_hier_sys(f'{top_dir}/top_DDI2_NTD.dat', name = 'tetramer2')

hier_sys, dof_sys, mols_sys = create_hier_sys(f'{top_dir}/topology.dat', name = 'system')

##############################
# Combined hierarchy
##############################

#p = IMP.Particle(mdl)
#hier_all = IMP.atom.Hierarchy.setup_particle(p)
#hier_all.add_child(hier_S1)
#hier_all.add_child(hier_S2)
#hier_all.set_name('System')

states = IMP.atom.get_by_type(hier_sys,IMP.atom.STATE_TYPE)
print('All states:', states)

##############################
# Connectivity
##############################
output_objects = [] 
sample_objects = []
rmf_restraints = []

crs = []

#add_connectivity(mols_S1, name='tetramer1')
#add_connectivity(mols_S2, name='tetramer2')
add_connectivity(mols_sys, name='system')
##############################
# Excluded Volume
##############################
#add_excluded_volume(mols_S1, name='tetramer1')
#add_excluded_volume(mols_S2, name='tetramer2')
add_excluded_volume(mols_sys, name='system')
##############################
# External barrier
##############################
#add_barrier(hier_S1, name = 'tetramer1')
#add_barrier(hier_S2, name = 'tetramer2')
add_barrier(hier_sys, name = 'system')
##############################
# XLs
##############################

#xls_tetramer = f'{top_dir}/data/ddi1_bifunctional.csv'
#add_amb_crosslinks(xls_tetramer, [hier_S1], dof_S1, weight = 5.0, length = 21., slope=0.01 ,name='tetramer1')
#
#xls_tetramer = f'{top_dir}/data/ddi2_bifunctional.csv'
#add_amb_crosslinks(xls_tetramer, [hier_S2], dof_S2, weight = 5.0, length = 21., slope=0.01 ,name='tetramer2')

xls_complete = f'{top_dir}/data/ddi_bifunctional.csv'
add_amb_crosslinks(xls_complete, [hier_sys], dof_sys, weight = 5.0, length = 21., slope=0.01 ,name='system')
##############################
# Distance Restraints from CSV
##############################
#add_distance_restraints_from_csv(xl_fil, hier_S1, min_distance=0.0, max_distance=19.0, resolution=1.0, kappa=1.0, weight=1.0, name='tetramer1')
#add_distance_restraints_from_csv(xl_fil, hier_S2, min_distance=0.0, max_distance=19.0, resolution=1.0, kappa=1.0, weight=1.0, name='tetramer2')
add_distance_restraints_from_csv(xl_fil, hier_sys, min_distance=0.0, max_distance=19.0, resolution=1.0, kappa=1.0, weight=1.0, name='system')
##############################
# Write everything
##############################  
    
output = IMP.pmi.output.Output()
output.init_rmf("ini_all.rmf3", [hier_sys])
output.write_rmf("ini_all.rmf3")
output.close_rmf("ini_all.rmf3")

##############################
# Shuffle
##############################    
#IMP.pmi.tools.shuffle_configuration(hier_S1, max_translation=100)
#IMP.pmi.tools.shuffle_configuration(hier_S2, max_translation=100)

IMP.pmi.tools.shuffle_configuration(hier_sys, max_translation=100)

movers_all = dof_sys.get_movers()

##############################
# MC
##############################    

rex=IMP.pmi.macros.ReplicaExchange(mdl,
                                   root_hier=hier_sys,                                        
                                   monte_carlo_sample_objects=movers_all,   
                                   replica_exchange_maximum_temperature=4.0,
                                   global_output_directory="output/",
                                   output_objects=output_objects,
                                   monte_carlo_steps=10,
                                   number_of_frames=1000,
                                   number_of_best_scoring_models=0)

rex.execute_macro()

exit()
