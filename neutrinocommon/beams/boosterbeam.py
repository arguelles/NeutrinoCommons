"""
Author  : C. Arguelles
Date    : 10/MAY/2011

This script models the booster beam.

Log:
- Modified on 23/FEB/2012 by C.Arguelles
    + Review the code to make it compatible with the neutrino
    commons library.
"""
import numpy as np
#my modules
import neutrinocommon.tools.generaltools as gt

# global variables
filepath = neutrinocommon.beams.__path__
global_datapath = filepath[0]+"/data/"

def LoadBoosterFlux(mbparam):
    """ Reads booster flux and stores it in memory.
    
    Ref. : A.A. Aguilar-Arevalo et al., "Measurement of Neutrino-Induced Charged Current-Charged Pion Production Cross Sections on Mineral Oil at Enu ~ 1 GeV"
    arXiv:1011.3572 [hep-ex]
    
    Data : http://www-boone.fnal.gov/for_physicists/data_release/ccpip/

    @type   mbparam   :   miniboone_config
    @param  mbparam   :   miniboone run configuration
    
    @rtype            :   array
    @return           :   raw flux data array
    """
    if mbparam.Beam.LoadFlux == [] or mbparam.Beam.LoadHornPolarity != mbparam.Beam.HornPolarity:
        if mbparam.Beam.HornPolarity == "positive":
            file = open(global_datapath + "pospolarity_fluxes.dat",'r')
        elif mbparam.Beam.HornPolarity == "negative":
            file = open(global_datapath + "negpolarity_fluxes.dat",'r')
        else:
            print "Wrong horn polarity."
            quit()
        h,dat = gt.hreadfilev2(file)
        file.close()
        mbparam.Beam.LoadFlux            = dat
        mbparam.Beam.LoadHornPolarity    = mbparam.Beam.HornPolarity
        return dat
    else :
        return mbparam.Beam.LoadFlux

def BoosterFlux(E,mbparam):
    """ Returns booster flux for a given horn configuration.
    
    Ref. : A.A. Aguilar-Arevalo et al., "Measurement of Neutrino-Induced Charged Current-Charged Pion Production Cross Sections on Mineral Oil at Enu ~ 1 GeV"
    arXiv:1011.3572 [hep-ex]

    @type   E         :   float
    @param  E         :   (anti)neutrino energy         [eV]
    @type   mbparam   :   miniboone_config
    @param  mbparam   :   miniboone run configuration
    
    @rtype            :   float
    @return           :   flux  [nu_X/cm^2/POT/eV]
    """
    flux_data = np.array(LoadBoosterFlux(mbparam))
    E_lo      = flux_data[:,0]*mbparam.GeV
    E_hi      = flux_data[:,1]*mbparam.GeV
    
    nu_mu     = flux_data[:,2]/(50.0*mbparam.MeV) # conv. scale to eV
    nu_mub    = flux_data[:,3]/(50.0*mbparam.MeV) # conv. scale to eV
    
    nu_e      = flux_data[:,4]/(50.0*mbparam.MeV) # conv. scale to eV
    nu_eb     = flux_data[:,5]/(50.0*mbparam.MeV) # conv. scale to eV
    
    for i,EE in enumerate(E_lo):
        if E >= E_lo[i] and E < E_hi[i]:
            return [nu_e[i],nu_eb[i],nu_mu[i],nu_mub[i]]
        else :
            pass
        
    return [0.0,0.0,0.0,0.0]
    
if __name__ == '__main__':
   # TESTING AREA #
   import minibooneclass as mbc
   mbcfg = mbc.MiniBooneExperiment()
   
   E = 200.0*mbcfg.MeV
   print BoosterFlux(E,mbcfg)
   
