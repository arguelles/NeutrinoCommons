""" 
Author  : C.A. Arguelles
Date    : 10/MAY/2011

This script constains the MiniBoone experiment information
and configuration. It is organized in several classes
for its easier setup.

Log :
- Modified on 23/ABR/2012 by C.Arguelles
    + Review the code to make it compatible with the neutrino
    commons library.  
"""
# python standard modules
import numpy as np
# my modules
import neutrinocommon.physconst.physicsconstants as PC

class MiniBooneExperiment(PC.PhysicsConstants):
    """
    MiniBoone Experiment Object
    """
        
    neutype = "neutrino"                    # make a neutrino or antineutrino analisis
    # Energy
    minenergy   = 0.2                       # [GeV]
    maxenergy   = 10.0                      # [GeV]
    stepenergy  = 0.2                       # [GeV]
        
    # OSCILATION CONFIGURATION
    class NeuOscConf():
        generations = True                  # If True : 3 generation formalism will be used. Else two generation formalism will be used.
        oscillation = True                  # If True : neutrino oscillations will be considered.
        matter 	    = True                  # If True : matter effects will be taken into account. Else vacuum oscillations will be used.)
        
    class Beam():
        LoadFlux            = []            # Loaded booster flux data      (aux)
        LoadHornPolarity    = ""            # Loaded booster horn polarity  (aux)
            
        HornPolarity        = "positive"    # Booster focusing horn polarity            
            
        POT_neutrino 	    = 6.46e20	    # Neutrino run Protons-on-target
        POT_antineutrino    = 8.5e20	    # Antineutrino run Protons-on-target
    
    class CrossSection():
        LoadCrossSection    = []            # Loaded cross section data     (aux)
        CalculationMode     = "binned"      # binned or interpolated CS
        
        ch      = {'CC_neutral_pion':0,'NC_neutral_pion':1,'CC_charged_pion':2,'CC_QE':3}
        chname = {0 : r"$\nu_\mu + n \rightarrow \mu^- + \pi_0 + p$",1 : r"$\nu_\mu + p(n) \rightarrow \nu_\mu + \pi_0 + p(n)$",2 : r"$\nu_\mu + p(n) \rightarrow \mu^- + \pi^+ + p(n)$",3 : r"$\nu_\mu + n \rightarrow \mu^- + p$"}
        
    class Detector():
        """ MiniBooNE detector specification.
        """
        Volume      = (4.0/3.0)*np.pi*(574.6)**3     #	[cm^3]
        Distance    = 541.0			     #	[m]
        
        class Material():
            Name            = 'CH_2'
            Density         = 0.845             #       [gr/cm^3]
            NucleonNum      = 14.0      	#	[nucleons/nuclei]
            AtomicWeight    = 14.0	        #	[gr/mol]
            
    class Events():
        def __init__(self):
            # electron neutrino events
            
            nu_e_evt_CC_QE           = 0.0
            nu_e_evt_CC_charged_pion = 0.0
            nu_e_evt_NC_neutral_pion = 0.0
            nu_e_evt_CC_neutral_pion = 0.0
                            
            # muon neutrino events                            
            nu_mu_evt_CC_QE           = 0.0
            nu_mu_evt_CC_charged_pion = 0.0
            nu_mu_evt_NC_neutral_pion = 0.0
            nu_mu_evt_CC_neutral_pion = 0.0
            
if __name__ == '__main__':
    pass
