from __future__ import print_function
import IMP
import IMP.core
import IMP.algebra
import IMP.atom
import IMP.container
import IMP.pmi.tools
import IMP.pmi.restraints
import numpy as np


class COMRestraint(IMP.pmi.restraints.RestraintBase):

    """ Positional restraint based on equivalence assignment """

    def __init__(self,
                 hier,
                 selection1,
                 selection2,
                 distance,
                 sigma_init = 6.0,
                 weight = 1.,
                 label = None):

        """Constructor:
        Generate  restraints of the absolute positions of
        proteins
        
        @param hier             The hierarchy 
        @param slope            Linear term to add to the restraint function to 
                                help when far away
        @param resolution       Beads resolution at which the restraint is applied
        
        """

        self.hier = hier
        self.selection1 = selection1
        self.selection2 = selection2
        self.sigma_init = sigma_init
        self.D = distance
        self.restraint_type = "D"
        
        model = self.hier.get_model()
        self.model = model
        
        rname = "COMRestraint"
        super(COMRestraint, self).__init__(
            model, name="COMRestraint", label=label, weight=weight)

        self.sigma_is_sampled = True
        self.sigma_dictionary={}

        self.rs = self.rs = IMP.RestraintSet(self.model, "COM_restraint")
        self.rs_sig = self._create_restraint_set("sigma")
        
        self._set_COM_restraint()
        
        #self.particles = particles

    def _get_particles_selection(self, selection):
        if isinstance(selection,tuple):
            prot =  selection[2]
            r1_1 = int(selection[0])
            r1_2 = int(selection[1])
            
            if len(selection)==4:
                ps1 = IMP.atom.Selection(self.hier,
                                         molecule=prot,
                                         residue_indexes=range(r1_1, r1_2),
                                         copy_index = selection[3]).get_selected_particles()
            else:
                ps1 = IMP.atom.Selection(self.hier,
                                         molecule=prot,
                                         residue_indexes=range(r1_1, r1_2)).get_selected_particles()
        
        else:
            ps1 = IMP.atom.Selection(self.hier,
                                     molecule=selection).get_selected_particles()

        return ps1

    def _set_COM_restraint(self):
                
        # Create nuisance particle
        sigma2name = 'sigma_pos'
        self.sigma = self._create_sigma(sigma2name, self.sigma_init)    
        # Add each assignment as an individual contribution
        
        self.ps1 = self._get_particles_selection(self.selection1)
        self.ps2 = self._get_particles_selection(self.selection2)

        

        cmr = IMP.isd.COMRestraint(self.model,
                                   self.D,
                                   self.sigma.get_particle_index(),
                                   self.restraint_type,
                                   False,
                                   'COMRestraint')

        cmr.add_contribution(self.ps1, '1')
        cmr.add_contribution(self.ps2, '2')

        self.rs.add_restraint(cmr)

        self.restraint_sets = [self.rs] + self.restraint_sets[1:]
        
    def _create_sigma(self, name, sigmainit):
        """ This is called internally. Creates a nuisance
        on the structural uncertainty """
        if name in self.sigma_dictionary:
            return self.sigma_dictionary[name][0]

        sigmaminnuis = 1.0
        sigmamaxnuis = 8.0
        sigmamin = 1.0
        sigmamax = 8.0
        sigmatrans = 0.2
        sigma = IMP.pmi.tools.SetupNuisance(self.model,
                                            sigmainit,
                                            sigmaminnuis,
                                            sigmamaxnuis,
                                            self.sigma_is_sampled).get_particle()
        self.sigma_dictionary[name] = (sigma,
                                       sigmatrans,
                                       self.sigma_is_sampled)
        
        self.rs_sig.add_restraint(IMP.isd.UniformPrior(self.model,
                                                       sigma,
                                                       100.0,
                                                       sigmamax,
                                                       sigmamin))
        
        self.rs_sig.add_restraint(IMP.isd.JeffreysRestraint(self.model, sigma))
        return sigma


    def get_particles_to_sample(self):
        ps = {}
        if self.sigma_is_sampled:
            for sigmaname in self.sigma_dictionary.keys():
                ps["Nuisances_COMRestraint_Sigma_" +
                   str(sigmaname) + self._label_suffix] = \
                    ([self.sigma_dictionary[sigmaname][0]],
                     self.sigma_dictionary[sigmaname][1])
        
        print('COM RESTRAINT> Number of nuisance particles:', len(ps))
        return ps
  

    '''
    def evaluate(self):

        self.ps1_centroid = IMP.core.get_centroid(IMP.core.XYZs(self.ps1))
        self.ps2_centroid = IMP.core.get_centroid(IMP.core.XYZs(self.ps2))

        print(self.ps1_centroid, type(self.ps1_centroid), self.ps2_centroid)
        
        dist = IMP.algebra.get_distance(self.ps1_centroid,
                                        self.ps2_centroid)

        print('dddd', dist)

        print(dist, self.D)
        if dist<=self.D:
            return 0.0
        else:
            print('-----', -np.log(np.exp(-(dist-self.D)*(dist-self.D)/8.0)))
            return -np.log(np.exp(-(dist-self.D)*(dist-self.D)/8.0))
    
    '''

    def get_output(self):

        output = super(COMRestraint, self).get_output()
        
        
        output["COMRestraint_Sigma" +
               self._label_suffix] = str(
                   self.sigma.get_scale())

        return output
