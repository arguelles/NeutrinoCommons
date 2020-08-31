""" 
Author  : C.A. Arguelles
Date    : 10/MAY/2011

This package contains functions to simulate the MiniBoone
experiment. This work was made using as a reference
the MiniBoone website information and associated papers.

Log :
- Modified on 23/ABR/2012 by C.Arguelles
    + Adapted the code to work with the neutrino
    commons package.
"""

# my modules
import neutrinocommon.astro.body as bd
import neutrinocommon.beams.boosterbeam as bb
import neutrinocommon.neu.xsection as xs
import neutrinocommon.neu.oscana as oa

#===============================================================================
# MINIBOONE DETECTOR EVENT ESTIMATION
#===============================================================================

def NucleonNumber(mbparam):
    """ Calculates the number of nucleons in MiniBooNE fiducial volume.
  
      @type   mbparam 	:	miniboone config. parameters
      @param  mbparam 	:	miniboone configuration information/parameters
    
      @rtype	        :	float
      @return	        :	number of nucleons
    """
    Det = mbparam.Detector
    DetMat = Det.Material
    
    Mass_eV = (Det.Volume*mbparam.cm**3)*(DetMat.Density*(mbparam.gr/mbparam.cm**3))
    
    return (Mass_eV/mbparam.gr)*mbparam.Na
    
def MiniBooNE_Efficiency(E,mbparam):
    """ Returns event efficiencies for different channels
  
      @type   E 	:	float
      @param  E 	:	energy  [eV]
      @type   mbparam 	:	miniboone config. parameters
      @param  mbparam 	:	miniboone configuration information/parameters
    
      @rtype	        :	float
      @return	        :	event efficiency in the channels
    """
    
    # Energy depedence is unknown
    # REFERENCES ???
    
    CC_QE           = 0.390
    CC_charged_pion = 0.457
    NC_neutral_pion = 1.000
    CC_neutral_pion = 0.064
    
    return [CC_neutral_pion,NC_neutral_pion,CC_charged_pion,CC_QE]
    
def FluxAtDetector(E,mbparam):
    """ Calculates the neutrino flux at the detector considering neutrino oscillations.
  
      @type   E 	:	float
      @param  E 	:	energy  [eV]
      @type   mbparam 	:	miniboone config. parameters
      @param  mbparam 	:	miniboone configuration information/parameters
    
      @rtype	        :	float
      @return	        :	flux [nu_X/cm^2/POT/eV]
    """    
    oscconf = mbparam.NeuOscConf
        
    if oscconf.oscillation :
        if oscconf.matter :
            body = bd.Earth()            
        else :
            body = bd.Vacuum()
        
        baseline    = mbparam.Detector.Distance*mbparam.meter
        track       = body.track(baseline,baseline,baseline)
        
        if oscconf.generations :
            mbparam.neutype == "neutrino"
            nu_e    = bb.BoosterFlux(E,mbparam)[0]*oa.NeuOsc3g(0,0,E,track,body,mbparam)+bb.BoosterFlux(E,mbparam)[2]*oa.NeuOsc3g(1,0,E,track,body,mbparam)
            nu_mu   = bb.BoosterFlux(E,mbparam)[2]*oa.NeuOsc3g(1,1,E,track,body,mbparam)+bb.BoosterFlux(E,mbparam)[0]*oa.NeuOsc3g(0,1,E,track,body,mbparam)
            mbparam.neutype == "antineutrino"
            nu_eb   = bb.BoosterFlux(E,mbparam)[1]*oa.NeuOsc3g(0,0,E,track,body,mbparam)+bb.BoosterFlux(E,mbparam)[3]*oa.NeuOsc3g(1,0,E,track,body,mbparam)
            nu_mub  = bb.BoosterFlux(E,mbparam)[3]*oa.NeuOsc3g(1,1,E,track,body,mbparam)+bb.BoosterFlux(E,mbparam)[1]*oa.NeuOsc3g(0,1,E,track,body,mbparam)
            return [nu_e,nu_eb,nu_mu,nu_mub]
        else :
            mbparam.neutype == "neutrino"
            nu_e    = bb.BoosterFlux(E,mbparam)[0]*oa.NeuOsc2g(0,0,E,track,body,mbparam)+bb.BoosterFlux(E,mbparam)[2]*oa.NeuOsc2g(1,0,E,track,body,mbparam)
            nu_mu   = bb.BoosterFlux(E,mbparam)[2]*oa.NeuOsc2g(1,1,E,track,body,mbparam)+bb.BoosterFlux(E,mbparam)[0]*oa.NeuOsc2g(0,1,E,track,body,mbparam)
            mbparam.neutype == "antineutrino"
            nu_eb   = bb.BoosterFlux(E,mbparam)[1]*oa.NeuOsc2g(0,0,E,track,body,mbparam)+bb.BoosterFlux(E,mbparam)[3]*oa.NeuOsc2g(1,0,E,track,body,mbparam)
            nu_mub  = bb.BoosterFlux(E,mbparam)[3]*oa.NeuOsc2g(1,1,E,track,body,mbparam)+bb.BoosterFlux(E,mbparam)[1]*oa.NeuOsc2g(0,1,E,track,body,mbparam)
            return [nu_e,nu_eb,nu_mu,nu_mub]
    else : 
        return bb.BoosterFlux(E,mbparam)
    
def NevtSpectraMiniBooNE(E,mbparam):
    """ Calculates the event number at the detector considering oscillations
     in every channel.
  
      @type   E 	:	float
      @param  E 	:	energy  [eV]
      @type   mbparam 	:	miniboone config. parameters
      @param  mbparam 	:	miniboone configuration information/parameters
    
      @rtype	        :	float
      @return	        :	event spectra [eV^-1]
    """
    if mbparam.Beam.HornPolarity == "positive":
        POT = mbparam.Beam.POT_neutrino
    elif mbparam.Beam.HornPolarity == "negative":
        POT = mbparam.Beam.POT_antineutrino
    else :
        print "Wrong horn polarization."
        quit()
        
    ch = mbparam.CrossSection.ch
    #{'CC_neutral_pion':0,'NC_neutral_pion':1,'CC_charged_pion':2,'CC_quasielastic':3}
        
    mbparam.neutype = "neutrino"
    
    evt_res = mbparam.Events()
    
    # electron neutrino events
    evt_res.nu_e_evt_CC_QE           = POT*NucleonNumber(mbparam)*xs.MiniBooNECS(E,mbparam)[ch['CC_QE']]*FluxAtDetector(E,mbparam)[0]
    evt_res.nu_e_evt_CC_charged_pion = POT*NucleonNumber(mbparam)*xs.MiniBooNECS(E,mbparam)[ch['CC_charged_pion']]*FluxAtDetector(E,mbparam)[0]
    evt_res.nu_e_evt_NC_neutral_pion = POT*NucleonNumber(mbparam)*xs.MiniBooNECS(E,mbparam)[ch['NC_neutral_pion']]*FluxAtDetector(E,mbparam)[0]
    evt_res.nu_e_evt_CC_neutral_pion = POT*NucleonNumber(mbparam)*xs.MiniBooNECS(E,mbparam)[ch['CC_neutral_pion']]*FluxAtDetector(E,mbparam)[0]
    
    # muon neutrino events    
    evt_res.nu_mu_evt_CC_QE           = POT*NucleonNumber(mbparam)*xs.MiniBooNECS(E,mbparam)[ch['CC_QE']]*FluxAtDetector(E,mbparam)[2]
    evt_res.nu_mu_evt_CC_charged_pion = POT*NucleonNumber(mbparam)*xs.MiniBooNECS(E,mbparam)[ch['CC_charged_pion']]*FluxAtDetector(E,mbparam)[2]#*MiniBooNE_Efficiency(E,mbparam)[ch['CC_charged_pion']]
    evt_res.nu_mu_evt_NC_neutral_pion = POT*NucleonNumber(mbparam)*xs.MiniBooNECS(E,mbparam)[ch['NC_neutral_pion']]*FluxAtDetector(E,mbparam)[2]#*MiniBooNE_Efficiency(E,mbparam)[ch['NC_neutral_pion']]
    evt_res.nu_mu_evt_CC_neutral_pion = POT*NucleonNumber(mbparam)*xs.MiniBooNECS(E,mbparam)[ch['CC_neutral_pion']]*FluxAtDetector(E,mbparam)[2]#*MiniBooNE_Efficiency(E,mbparam)[ch['CC_neutral_pion']]
    
    #e_evt   = [nu_e_evt_CC_QE,nu_e_evt_CC_charged_pion,nu_e_evt_NC_neutral_pion,nu_e_evt_CC_charged_pion]
    #mu_evt  = [nu_mu_evt_CC_QE,nu_mu_evt_CC_charged_pion,nu_mu_evt_NC_neutral_pion,nu_mu_evt_CC_charged_pion]
    
    #neutrino_evt = [e_evt,mu_evt]
    
    #mbparam.neutype = "antineutrino"
    #
    #evt_nub_CC_QE           = POT*NucleonNumber(mbparam)*(xs.MiniBooNECS(E,mbparam)[ch['CC_QE']]*FluxAtDetector(E,mbparam)[1]+xs.MiniBooNECS(E,mbparam)[ch['CC_QE']]*FluxAtDetector(E,mbparam)[3])*MiniBooNE_Efficiency(E,mbparam)[ch['CC_QE']]
    #evt_nub_CC_charged_pion = POT*NucleonNumber(mbparam)*(xs.MiniBooNECS(E,mbparam)[ch['CC_charged_pion']]*FluxAtDetector(E,mbparam)[3])*MiniBooNE_Efficiency(E,mbparam)[ch['CC_charged_pion']]
    #evt_nub_NC_neutral_pion = POT*NucleonNumber(mbparam)*(xs.MiniBooNECS(E,mbparam)[ch['NC_neutral_pion']]*FluxAtDetector(E,mbparam)[3])*MiniBooNE_Efficiency(E,mbparam)[ch['NC_neutral_pion']]
    #evt_nub_CC_neutral_pion = POT*NucleonNumber(mbparam)*(xs.MiniBooNECS(E,mbparam)[ch['CC_neutral_pion']]*FluxAtDetector(E,mbparam)[3])*MiniBooNE_Efficiency(E,mbparam)[ch['CC_neutral_pion']]
    
    #antineutrino_evt    = [evt_nub_CC_neutral_pion,evt_nub_NC_neutral_pion,evt_nub_CC_charged_pion,evt_nub_CC_QE]
    antineutrino_evt = []
        
    #return [neutrino_evt,antineutrino_evt]
    return evt_res
    
if __name__ == '__main__':
    # TESTING AREA #
    import minibooneclass as mbc
    mbparam = mbc.MiniBooneExperiment()
    
    E = 1.0*mbparam.GeV
    
    print NevtSpectraMiniBooNE(E,mbparam)
    
    
    
    
