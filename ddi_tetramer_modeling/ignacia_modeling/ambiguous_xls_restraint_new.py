"""@namespace IMP.pmi.restraints.pemap
Restraints for handling pemap data.
"""

#!/usr/bin/env python
import IMP
import IMP.core
import IMP.algebra
import IMP.atom
import IMP.isd
import IMP.container
import IMP.pmi.tools
import IMP.pmi.output
import IMP.pmi.restraints
from math import log
from collections import defaultdict
import itertools
import operator
import os
import math
import pandas as pd
import numpy as np
    
class XLs_amb_Restraint(IMP.pmi.restraints.RestraintBase):

    _include_in_rmf = True
    def __init__(self,
                 hiers,
                 xls_file,
                 sigma_init = 4.0,
                 length = 30.,
                 slope = 0.0,
                 weight=1.0,
                 label = None):

        """
        Constructor:
        @param [hiers]           Hierarchies
        @param xls_file          File containing p1,r1,p2,r2,id
        @param sigma_init        Initial value for noise parameter
        @param weight            Weight of the restraint
        @param label             Extra text to label the restraint so that it is
                                 searchable in the output  
        """

        self.root_hiers = hiers
        if isinstance(hiers, list): 
            self.m = hiers[0].get_model()
        else:
            self.m = hiers.get_model()
            self.root_hiers = [hiers]
            
        rname = "CrossLinkingAmbiguosRestraint"
        super(XLs_amb_Restraint, self).__init__(
            self.m, name=rname, label=label, weight=weight)

        self.xls_file = xls_file

        #self.length = length
        #self.slope = slope
        #self.weight = weight

       
        # Restraint information
        self.rs.set_name(self.rs.get_name() + "_Data")
        self.rspsi = self._create_restraint_set("PriorPsi")
        self.rssigma = self._create_restraint_set("PriorSigma")
        self.rslin = self._create_restraint_set("Linear")
        # Add custom metadata (will be saved in RMF output)
        self.rs.length = length
        self.rs.slope = slope
        

        # dummy linear restraint used for Chimera display
        self.linear = IMP.core.Linear(0, 0.0)
        self.linear.set_slope(0.0)
        self.dps2 = IMP.core.DistancePairScore(self.linear)
        
        print('-------',  self.rs.length)
    
        # Nuisance particle parameters
        self.sigma_init = sigma_init
        self.psi_is_sampled = True
        self.sigma_is_sampled = False

        self.psi_dictionary={}
        self.sigma_dictionary={}

        # Create unique sigma parameter
        self.sigma = self._create_sigma('sigma', self.sigma_init)        
        
        # Read file
        self._read_xls_file()

        # Create restraint
        self._create_restraints()

    def _create_restraints(self):
        # Setup restraint
        
        self.xl_list = []
        restraints = []
        print('all_xls', len(self.all_xls))
        print('all_xls', self.all_xls)
    

        for id, sys_dict in self.all_xls.items():
            c= 0
            dr = IMP.isd.CrossLinkMSAmbRestraint(self.model,
                                                 self.rs.length,
                                                 self.rs.slope)
            xl_label = id
            for sys, xl in sys_dict.items():
                psi = xl[0][2].get_particle_index()
                pairs = [[p[0].get_index(),p[1].get_index()] for p in xl]
                dr.add_contribution(sys,
                                    pairs,
                                    self.sigma.get_particle_index(),
                                    psi)

                for pi1, pi2 in pairs:
                    pr = IMP.core.PairRestraint(self.model, self.dps2,
                                                (pi1, pi2))
                    pr.set_name(str(xl_label))
                    self.rslin.add_restraint(pr)
                    restraints.append(dr)
                
        if len(restraints) == 0:
            raise SystemError(
                "CrossLinkingMassSpectrometryRestraint: no cross-link "
                "was constructed")
        self.xl_restraints = restraints
        lw = IMP.isd.LogWrapper(restraints, 1.0)
        self.rs.add_restraint(lw)
    
    def _read_xls_file(self):
        d = pd.read_csv(self.xls_file)
        gd = d.groupby('Iid')

        self.labels = {}
        
        # Pairs of crosslinked beads for each system
        self.all_xls = {}
        
        for i, (name, group) in enumerate(gd):
            id = group.iloc[0]["Iid"]
            self.labels[i]= id
            g_sel = {j:[] for j in range(len(self.root_hiers))}

            # Create PSI parameters
            if 'Score' in group.columns:
                psiname = f'PSI_{group.iloc[0]["Score"]}'
                psi = self.create_psi(psiname)
            else:
                psiname = f'PSI'
                psi = self.create_psi(psiname)
            
            # Find beads  
            for m, row in group.iterrows():
                for n, h in enumerate(self.root_hiers):
                    # Select particles
                    s1,s2 = self._select_particles(h, row)
                    
                    if len(s1)>0 and len(s2)>0:
                        for pair in itertools.product(s1,s2):
                            if pair[0] != pair[1]:
                                g_sel[n].append([pair[0],pair[1],psi])
                            
            # Linkages found in representation
            g_sel = {k:v for k,v in g_sel.items() if len(v)>0}

            if len(g_sel)>0:
                self.all_xls[name] = g_sel

        print('------', self.all_xls)
                
    def _select_particles(self, hier, row):
        if 'Copy A' in row.index and 'Copy B' in row.index:
            copy_index_A = row['Copy A']
            copy_index_B = row['Copy B']
            s1 = IMP.atom.Selection(hier,
                                    molecule = row['Prot A'],
                                    residue_index = int(row['Residue A']),
                                    copy_index=copy_index_A).get_selected_particles()
            s2 = IMP.atom.Selection(hier,
                                    molecule = row['Prot B'], 
                                    residue_index = int(row['Residue B']),
                                    copy_index=copy_index_B).get_selected_particles()
            

        else:
            s1 = IMP.atom.Selection(hier,
                                    molecule = row['Prot A'],
                                    residue_index = int(row['Residue A'])).get_selected_particles()
            s2 = IMP.atom.Selection(hier,
                                    molecule = row['Prot B'], 
                                    residue_index = int(row['Residue B'])).get_selected_particles()

        return s1, s2


    def create_xi(self, name):
        """ This is called internally. Creates a nuisance
        on the data uncertainty """
        if name in self.xi_dictionary:
            return self.xi_dictionary[name][0]

        xiinit = 0.5
        ximinnuis = 0.0000001
        ximaxnuis = 0.9999999
        ximin = 0.01
        ximax = 0.99
        xitrans = 0.1
        xi = IMP.pmi.tools.SetupNuisance(
            self.model, xiinit, ximinnuis, ximaxnuis,
            self.xi_is_sampled, name=name).get_particle()
        self.xi_dictionary[name] = (xi, xitrans, self.xi_is_sampled)

        self.rsxi.add_restraint(
            IMP.isd.UniformPrior(self.model,
                                 xi,
                                 10000.0,
                                 ximax,
                                 ximin))

        self.rsxi.add_restraint(IMP.isd.JeffreysRestraint(self.model, xi))
        return xi       
    
    def create_psi(self, name):
        """ This is called internally. Creates a nuisance
        on the data uncertainty """
        if name in self.psi_dictionary:
            return self.psi_dictionary[name][0]

        psiinit = 0.25
        psiminnuis = 0.0000001
        psimaxnuis = 0.4999999
        psimin = 0.01
        psimax = 0.49
        psitrans = 0.1
        psi = IMP.pmi.tools.SetupNuisance(
            self.model, psiinit, psiminnuis, psimaxnuis,
            self.psi_is_sampled, name=name).get_particle()
        self.psi_dictionary[name] = (psi, psitrans, self.psi_is_sampled)

        self.rspsi.add_restraint(
            IMP.isd.UniformPrior(self.model, psi, 1000000000.0,
                                 psimax, psimin))

        self.rspsi.add_restraint(IMP.isd.JeffreysRestraint(self.model, psi))
        return psi
            
    def _create_sigma(self, name, sigmainit):
        """ This is called internally. Creates a nuisance
        on the structural uncertainty """
        if name in self.sigma_dictionary:
            return self.sigma_dictionary[name][0]

        #m = self.root_hier.get_model()
        
        sigmainit = 2.0
        sigmaminnuis = 0.0000001
        sigmamaxnuis = 1000.0
        sigmamin = 0.01
        sigmamax = 100.0
        sigmatrans = 0.5
        sigma = IMP.pmi.tools.SetupNuisance(
            self.model, sigmainit, sigmaminnuis, sigmamaxnuis,
            self.sigma_is_sampled, name=name).get_particle()
        self.sigma_dictionary[name] = (
            sigma,
            sigmatrans,
            self.sigma_is_sampled)
        self.rssigma.add_restraint(
            IMP.isd.UniformPrior(
                self.model,
                sigma,
                1000000000.0,
                sigmamax,
                sigmamin))
        return sigma
    
    def get_particles_to_sample(self):
        """ Get the particles to be sampled by the IMP.pmi.sampler object """
        ps = {}
        if self.sigma_is_sampled:
            for sigmaname in self.sigma_dictionary:
                ps["Nuisances_CrossLinkinMSAmbRestraint_Sigma_" +
                   str(sigmaname) + self._label_suffix] =\
                    ([self.sigma_dictionary[sigmaname][0]],
                     self.sigma_dictionary[sigmaname][1])

        if self.psi_is_sampled:
            for psiname in self.psi_dictionary:
                ps["Nuisances_CrossLinkingMSAmbRestraint_Psi_" +
                   str(psiname) + self._label_suffix] = \
                   ([self.psi_dictionary[psiname][0]],
                    self.psi_dictionary[psiname][1])
                
        print('AMB XLs RESTRAINT> Number of nuisance particles:', len(ps))
        
        return ps

    def get_movers(self):
        """ Get all need data to construct a mover in IMP.pmi.dof class"""
        movers = []
        if self.sigma_is_sampled:
            for sigmaname in self.sigma_dictionary:
                mover_name = \
                    "Nuisances_CrossLinkingMSAmbRestraint_Sigma_" \
                    + str(sigmaname) + "_" + self.label
                particle = self.sigma_dictionary[sigmaname][0]
                maxstep = (self.sigma_dictionary[sigmaname][1])
                mv = IMP.core.NormalMover(
                    [particle], IMP.FloatKeys([IMP.FloatKey("nuisance")]),
                    maxstep)
                mv.set_name(mover_name)
                movers.append(mv)

        if self.psi_is_sampled:
            for psiname in self.psi_dictionary:
                mover_name = \
                    "Nuisances_CrossLinkingMSAmbRestraint_Psi_" + \
                    str(psiname) + "_" + self.label
                particle = self.psi_dictionary[psiname][0]
                maxstep = (self.psi_dictionary[psiname][1])
                mv = IMP.core.NormalMover(
                    [particle], IMP.FloatKeys([IMP.FloatKey("nuisance")]),
                    maxstep)
                mv.set_name(mover_name)
                movers.append(mv)

        return movers

    def get_output(self):
        """Get the output of the restraint to be used by the IMP.pmi.output
           object"""
        output = super(XLs_amb_Restraint,
                       self).get_output()

        for psiname in self.psi_dictionary:
            output["CrossLinkingMSAmbRestraint_Psi_" +
                   str(psiname) + self._label_suffix] = str(
                        self.psi_dictionary[psiname][0].get_scale())

        for sigmaname in self.sigma_dictionary:
            output["CrossLinkingMSAmbRestraint_Sigma_" +
                   str(sigmaname) + self._label_suffix] = str(
                    self.sigma_dictionary[sigmaname][0].get_scale())

        return output

    
    def get_restraint_for_rmf(self):
        """ get the dummy restraints to be displayed in the rmf file """
        return self.rs, self.rssigma, self.rspsi
    
    
