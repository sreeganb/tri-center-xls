#####################################
# Modeling of Pom34-Pom152 TMDs
# using two-fold symmetry
#
# Ignacia Echeverria - Sali Lab
######################################


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
import IMP.pmi.restraints.em 
import IMP.pmi.restraints.basic
import IMP.pmi.io.crosslink
import IMP.pmi.restraints.npc
import IMP.pmi.mmcif
import IMP.bayesianem
import IMP.bayesianem.restraint
import math
import time
import numpy as np
import os
import sys
import ihm.cross_linkers

from npc import *
sys.path.append('../utils/')
import make_archive


top_dir = '../'

# Create System and State
mdl = IMP.Model()
s = IMP.pmi.topology.System(mdl)
st = s.create_state()

mols = []
domains = {'Pom152':[[105,130], [144,167], [176, 192], [200, 212]],
           'Pom34':[[44,86],[89,110],[122,150],[222,237]]}


pdb = {'Pom152':'data/AF-P39685-F1-model_v4.pdb',
       'Pom34':'data/AF-Q12445-F1-model_v4.pdb'}

chain_id = 'A'
n_densities = 10
colors = ['gray','red' ]

for k, (prot, doms) in enumerate(domains.items()):
    seqs = IMP.pmi.topology.Sequences(os.path.join(top_dir,f'data/{prot}.fasta'))
    mol = st.create_molecule(prot,sequence=seqs[prot],chain_id='A')
    pdb_file = os.path.join(top_dir,pdb[prot])
    struct = []
    for i, d in enumerate(doms):
        a1 = mol.add_structure(pdb_file,
                           chain_id=chain_id,
                           res_range = d)

        # Add structured part representation and then build
        mol.add_representation(a1,
                               density_residues_per_component=n_densities,
                               density_voxel_size=3.0,
                               resolutions=[1],
                               density_prefix = os.path.join(f'{top_dir}/data/em_data',f'{prot}_{i}'),
                               color = colors[k])
        struct += a1
    mol.add_representation(mol[0:250]-struct,resolutions=[2], color = colors[k])
    mols.append(mol)
        
# Clone everything
clones = []
chains='DE'
for i,mol in enumerate(mols):
    clone = mol.create_clone(chains[i])
    clones.append(clone)

mols = mols + clones

############################################
# Generate mmcif file
############################################
    
if '--mmcif' in sys.argv:
    # Record the modeling protocol to an mmCIF file
    po = IMP.pmi.mmcif.ProtocolOutput()
    po.system.title = ('Integrative structure of the Pom34-Pom152 TMDs')
    s.add_protocol_output(po)

############################################
# Build hierarchy
############################################
hier = s.build()

mdl.update() # propagates coordinates

out = IMP.pmi.output.Output()
out.init_rmf("sym_ini.rmf3", [hier])
out.write_rmf("sym_ini.rmf3")
out.close_rmf("sym_ini.rmf3")

############################################
# DOFs
############################################             
dof = IMP.pmi.dof.DegreesOfFreedom(mdl)

for i, mol in enumerate(mols):
    dof.create_flexible_beads(mol.get_non_atomic_residues())
    copy_index =  IMP.atom.Copy(mol.get_hierarchy()).get_copy_index()
    print(mol.get_name(), copy_index)
    for dom in domains[mol.get_name()]:
        sel = IMP.atom.Selection(hier,
                                 molecule=mol.get_name(),
                                 residue_indexes=range(dom[0],dom[1]+1,1),
                                 copy_index=copy_index,
                                 resolution=IMP.atom.ALL_RESOLUTIONS).get_selected_particles()
        sel_densities = IMP.atom.Selection(hier,
                                           molecule=mol.get_name(),
                                           residue_indexes=range(dom[0],dom[1]+1,1),
                                           copy_index=copy_index,
                                           representation_type=IMP.atom.DENSITIES).get_selected_particles()

        dof.create_rigid_body(sel+sel_densities)
   
print(dof.get_movers())

############################################
# Symmetry
############################################

center = IMP.algebra.Vector3D([0,0,0])

# Axial symmetry
rot = IMP.algebra.get_rotation_about_axis([0,1,0],math.pi)
transform = IMP.algebra.get_rotation_about_point(center,rot)
dof.constrain_symmetry(mols[0],mols[2],transform)
dof.constrain_symmetry(mols[1],mols[3],transform)

############################################
# Connectivity restraint
############################################
output_objects = []

for n, mol in enumerate(mols):
    cr = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(mol)
    cr.set_label(mol.get_name()+'.'+str(n))
    cr.add_to_model()
    output_objects.append(cr)

############################################
# Membrane binding
############################################
tor_th      = 45.0
tor_th_ALPS = 12.0
tor_R       = 390.0 + 180.0
tor_r       = 150.0 - tor_th/2.0
tor_r_ALPS  = 150.0 - tor_th_ALPS/2.0
msl_sigma   = 2.0
msl_weight  = 1.0

############################################
# Membrane binding
############################################

transmembrane_sel = [(106,133,'Pom152',0),
                     (144,167,'Pom152',0),
                     (176,194,'Pom152',0),
                     (57,87,'Pom34',0),
                     (122,150,'Pom34',0)]

for sel in transmembrane_sel:
    msl = IMP.pmi.restraints.npc.MembraneSurfaceLocationRestraint(hier=hier,
                                                                  protein=sel,
                                                                  tor_R=tor_R,
                                                                  tor_r=tor_r,
                                                                  tor_th=tor_th,
                                                                  sigma=msl_sigma,
                                                                  resolution = 1)
    msl.set_weight(msl_weight)
    msl.add_to_model()
    #msl.set_label('%s.%s'%(sel[2],sel[0]))
    output_objects.append(msl)
    print('Membrane binding restraint:', msl.evaluate())

perinuclear_sel = [(134,143,'Pom152',0),
                   (195,210,'Pom152',0),
                   (88,89,'Pom34',0),
                   (111,118,'Pom34',0)]
                      
'''    
for sel in perinuclear_sel:
    print('Applying membrane localization restraint:', sel)
    msl_1 = PerinuclearVolumeLocationRestraint(hier,
                                                                      protein=sel,
                                                                      tor_R=tor_R,
                                                                      tor_r=tor_r,
                                                                      tor_th=tor_th,
                                                                      sigma=msl_sigma,
                                                                      resolution = 10,
                                                                      label=f'{sel[2]}.{sel[0]}')
    msl_1.set_weight(msl_weight)
    msl_1.add_to_model()
    output_objects.append(msl_1)
    print('Membrane binding restraint ready', msl_1.get_output())

poreside_sel = [(1,43,'Pom34',0),
                (158,299,'Pom34',0),
                (1,104,'Pom152',0),
                (167,175,'Pom152',0)]
    
    
for sel in poreside_sel:
    print('Applying membrane localization restraint:', sel)
    msl_1 = PoreSideVolumeLocationRestraint(hier,
                                                                   protein=sel,
                                                                   tor_R=tor_R,
                                                                   tor_r=tor_r,
                                                                   tor_th=tor_th,
                                                                   sigma=msl_sigma,
                                                                   resolution = 10,
                                                                   label=f'{sel[2]}.{sel[0]}')
    msl_1.set_weight(msl_weight)
    msl_1.add_to_model()
    output_objects.append(msl_1)
    print('Membrane binding restraint ready', msl_1.get_output())
'''    
############################################
# EM
############################################     

densities = IMP.atom.Selection(hier,
                               representation_type=IMP.atom.DENSITIES).get_selected_particles()
    
print('Number of densities', len(densities))
gem = IMP.bayesianem.restraint.GaussianEMRestraintWrapper(densities,
                                                          os.path.join(top_dir,'data/em_data/run27h-c2-50-300bf_locres_filt-mask2-2-norm-100-zone16-TMDs-dust_dsfact1_ng60_oriented.txt'),
                                                          scale_target_to_mass=True)
gem.set_label("EM_membrane")
gem.add_to_model()
gem.set_weight(20.0)
output_objects.append(gem)


t0 = gem.evaluate()

############################
# Excluded Volume
############################

# Create excluded volume for all particles
evr1 = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(included_objects=mols)
evr1.add_to_model()
evr1.set_weight(1.0)
#evr1.set_label('intra')
output_objects.append(evr1)

##########################
# Inter-molecular XLs
##########################
resi_XLs = [62, 301, 351]
resi_XLs = [(62,'Pom152'), (44,'Pom34')]
w_intra_xls = 5.
include_intra_XLs = True
if include_intra_XLs:
    for (r, prot) in resi_XLs:
        dist_min = 10.0
        dist_max = 30.0
        ixl = IMP.pmi.restraints.basic.DistanceRestraint(root_hier = hier,
                                                         tuple_selection1=(r,r,prot,0),
                                                         tuple_selection2=(r,r,prot,1),
                                                         distancemin=dist_min,
                                                         distancemax=dist_max,
                                                         label=f'XLs_inter_{r}')
        ixl.set_weight(w_intra_xls)
        ixl.add_to_model()
        output_objects.append(ixl)
        print('Intra molecular XLs:', ixl.get_output())

###########################
# Chemical crosslinks
###########################
# INITIALIZE DB
rmf_restraints = []
include_XLs = True
w_xls = 5.
if include_XLs:
    cldbkc=IMP.pmi.io.crosslink.CrossLinkDataBaseKeywordsConverter()
    cldbkc.set_protein1_key("Protein1")
    cldbkc.set_protein2_key("Protein2")
    cldbkc.set_residue1_key("Residue1")
    cldbkc.set_residue2_key("Residue2")
    
    # XLs RESTRAINT
    cldb=IMP.pmi.io.crosslink.CrossLinkDataBase(cldbkc)
    cldb.create_set_from_file(os.path.join(top_dir,"data/XLs_all_2020.csv"))

    xl1 = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(root_hier=hier,
                                                                                database=cldb,
                                                                                resolution=1,
                                                                                length=21.0,
                                                                                slope=0.01,
                                                                                linker=ihm.cross_linkers.dss)
    xl1.add_to_model()
    xl1.set_weight(w_xls)
    rmf_restraints.append(xl1)
    output_objects.append(xl1)
        
###########################
# Randomize configurations
###########################

sel = IMP.atom.Selection(hier).get_selected_particles()

IMP.pmi.tools.shuffle_configuration(sel,
                                    bounding_box=((-150, -540, -20), (150, -440, 20)),
                                    avoidcollision_rb=False)

mdl.update()

############################
# Sampling
############################

nframes = 100000
if '--mmcif' in sys.argv:
    nframes = 5

mc1 = IMP.pmi.macros.ReplicaExchange(mdl,
                                      root_hier=hier,                              
                                      monte_carlo_sample_objects=dof.get_movers(),  
                                      global_output_directory="output/",
                                      output_objects=output_objects,
                                      replica_exchange_maximum_temperature=4.0,
                                      monte_carlo_steps=20,
                                      number_of_frames=nframes,
                                      number_of_best_scoring_models=0)

mc1.execute_macro()
rex1 = mc1.get_replica_exchange_object()

##############################
# Generate mmcif
##############################  
if '--mmcif' in sys.argv:
    import ihm.cross_linkers
    import ihm.dumper
    import ihm.format
    import ihm.location
    import ihm.representation
    import ihm.startmodel
    import ihm.dataset
    import ihm.protocol
    import ihm.analysis
    import ihm.model
    import ihm.restraint
    import ihm.geometry
    
    fname = os.path.join(top_dir,"data/XLs_all_2020.csv")

    # Add publication
    po.system.title = "Implications of a multiscale structure of the yeast Nuclear Pore Complex"
    po.system.citations.append(ihm.Citation(
        pmid='NA',
        title="Implications of a multiscale structure of the yeast Nuclear Pore Complex",
        journal="Molecular Cell", 
        year=2023,
        volume='NA',
        doi='NA',
        page_range='NA',
        authors=['Akey CA', 'Echeverria I', 'Ouch C',
                 'Nudelman I', 'Shi Y', 'Wang J', 'Weiss TM', 'Shi Y',
                 'Chait BT', 'Sali A','Fernandez-Martinez J', 'Rout MP']))
    
    s = po.system
    print("restraint datasets:", [r.dataset for r in s.restraints])
    # Datasets for XL-MS restraint
    l = ihm.location.Repository(url='https://zenodo.org/record/5721514',
                                doi='10.5281/zenodo.5721514',
                                details='All raw mass spectrometry files and peaklists used in the study')
    for r in s.restraints:
        print('----', r)
        if isinstance(r, ihm.restraint.CrossLinkRestraint):
            r.linker = ihm.cross_linkers.dsso
            r.dataset.location.details = l
            print("XL-MS dataset at:", r.dataset.location.path)
            print("Details:", r.dataset.location.details)
    # Correct number of output models to account for multiple runs
    protocol = s.orphan_protocols[-1]
    protocol.steps[-1].num_models_end = 6400000

    # Get last protocol in the file
    protocol = po.system.orphan_protocols[-1]
    # State that we filtered the 200000 frames down to one cluster of
    # 9999 models:
    analysis = ihm.analysis.Analysis()
    protocol.analyses.append(analysis)
    analysis.steps.append(ihm.analysis.ClusterStep(
                            feature='RMSD',
                           num_models_begin=6400000,
                            num_models_end=7089))

    # Create an ensemble for the cluster
    #e = po._add_simple_ensemble(analysis.steps[-1],
    #                            name="Cluster 0", num_models=9999,
    #                            drmsd=8.3, num_models_deposited=1,
    #                            localization_densities={}, ensemble_file=None)

    mg = ihm.model.ModelGroup(name="Cluster 0")
    po.system.state_groups[-1][-1].append(mg)
    e = ihm.model.Ensemble(model_group=mg,
                       num_models=11,
                       post_process=analysis.steps[-1],
                       precision=9.52,
                       name="Cluster 0")
    po.system.ensembles.append(e)
    
    # Add the model from RMF - centroid
    rh = RMF.open_rmf_file_read_only('../results/clustering/cluster.0/cluster_center_model.rmf3')
    IMP.rmf.link_hierarchies(rh, [hier])
    IMP.rmf.load_frame(rh, RMF.FrameID(0))
    mdl.update()
    del rh
    model = po.add_model(e.model_group)

    # Add ensemble members
    models = []
    rmf_ensemble = '../results/cluster0_random_sel.rmf3'
    rh = RMF.open_rmf_file_read_only(rmf_ensemble)
    IMP.rmf.link_hierarchies(rh, [hier])
    for frame_number in range(rh.get_number_of_frames()):        
        IMP.rmf.load_frame(rh, RMF.FrameID(frame_number))
        mdl.update()
        models.append(hier)
        model = po.add_model(e.model_group)

    #model_group = ihm.model.ModelGroup(models, name="Cluster 0")
    #state = ihm.model.State([model_group])
    #s.state_groups.append(ihm.model.StateGroup([state]))
    
    #del rh
    
    # Add localization densities
    # Look up the ihm.AsymUnit corresponding to a PMI component name
    for asym in po.asym_units:
        name = asym.split('.')[0]
        fname = f'../results/clustering/cluster.0/LPD_{name}_0.mrc'
        print('fname', fname)
        loc = ihm.location.OutputFileLocation(fname)
        den = ihm.model.LocalizationDensity(file=loc, asym_unit=po.asym_units[asym])
        # Add to ensemble
        e.densities.append(den)
    
    # Add uniprot of proteins

    lpep = ihm.LPeptideAlphabet()
    d = 'Segments used for modeling'
    
    # name : (uniprot id, mutations, [[db_begin, db_end, entity_begin, entity_end]]
    Uniprot={'Pom34.0': ('Q12445',[1,299,1,250]),
             'Pom152.0': ('P39685',[1,1337,1,250])}
    
    for prot, (entry, limits) in Uniprot.items():
        print(prot, entry, limits)
        ref = ihm.reference.UniProtSequence.from_accession(entry)
        ref.alignments.append(ihm.reference.Alignment(
            db_begin=limits[0], db_end=limits[1], entity_begin=limits[2], entity_end=limits[3]))
            
        po.asym_units[prot].entity.references.append(ref)

    # Point to the EM density map
    lem = ihm.location.EMDBLocation('EMD-41117',
                                 details='Pom34-152 membrane attachment site in the yeast NPC')
    gem.dataset.add_primary(ihm.dataset.EMDensityDataset(location=lem))

     # Dataset for EM restraint
    em, = [r for r in s.restraints
           if isinstance(r, ihm.restraint.EM3DRestraint)]
    d = em.dataset
    print("GMM file at", d.location.path)
    print("is derived from EMDB entry", d.parents[0].location.access_code)

    # Replace local links with DOIs
    repos = []
    for subdir, zipname in make_archive.ARCHIVES.items():
        print('subdir', subdir)
        repos.append(ihm.location.Repository(
            doi="10.5281/zenodo.8226857", root="../%s" % subdir,
            url="https://zenodo.org/record/8226857/files/%s.zip" % zipname,
            top_directory=os.path.basename(subdir)))

    po.system.update_locations_in_repositories(repos)
    
    po.finalize()
    with open('IM_TMDs_dimer.cif', 'w') as fh:
        ihm.dumper.write(fh, [po.system])
    


exit()
