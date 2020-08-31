""" 
Author  : C.A. Arguelles
Date    : 10/MAY/2011

This package contains functions to simulate the MINOS 
experiment. This work was made using as a reference
the MINOS website information and associated papers.

References :

- Interpretation of MINOS data in terms of non-standard
neutrino interactions
J. Kopp, P. Machado, S. Parke
arXiv : 1009.0014
Appendix A

- Neutrino and antineutrino inclusive charged-current
cross section measurements with the MINOS near detector
The MINOS collaboration
arXiv : 0910.2201

Log :
- Modified on 23/ABR/2012 by C.Arguelles
    + Added referenced.
    + Adapted the code to work with the neutrino
    commons package.
"""

import numpy as np
import scipy.integrate as integrate
import scipy.interpolate as interpolate
# my modules
import neutrinocommon
import neutrinocommon.beams.numibeam as nb
import neutrinocommon.neu.xsections as xs
import neutrinocommon.neu.oscana as oa
import neutrinocommon.neu.neuosc as no
import neutrinocommon.astro.body as bd
import neutrinocommon.tools.generaltools as gt

# global variables
filepath = neutrinocommon.exp.minos.__path__
global_datapath = filepath[0]+"/data/"

#===============================================================================
# MINOS DETECTOR EVENT ESTIMATION
# Work based on : arXiv : 0910.2201
#===============================================================================

def TargetNumber(mparam):
    """ Calculates the number of nucleons in the MINOS far detector.
  
      @type   mparam 	:	minos config. parameters
      @param  mparam 	:	minos configuration information/parameters
    
      @rtype	        :	float
      @return	        :	number of nucleons
    """
    FarDet = mparam.FarDet
    FarDetMat = FarDet.Material
    return FarDetMat.NucleonNum*(FarDet.mass*mparam.kg)/(FarDetMat.AtomicWeight*mparam.gr)*mparam.Na
        
def MINOSEfficiency(E,mparam):
    """ Returns CC-event efficiencies for MINOS far detector.
  
      @type   E 	:	float
      @param  E 	:	energy  [eV]
      @type   mparam 	:	minos config. parameters
      @param  mparam 	:	minos configuration information/parameters
    
      @rtype	        :	float
      @return	        :	CC-event efficiency
    """
    file  = open(global_datapath + "MINOS_CC_efficiency.dat",'r')
    h,dat = gt.hreadfilev2(file)
    
    Enu   = np.array(dat)[:,0]
    Eff   = np.array(dat)[:,1]
    if E/mparam.GeV <= Enu[-1]:
        eff_interpolate = interpolate.interp1d(Enu,Eff)
        return eff_interpolate(E/mparam.GeV)
    else :
        return Eff[-1]
    
def FluxAtFarDetector(E,mparam):    
    """ Calculates flux spectra at MINOS far detector.
  
      @type   E 	:	float
      @param  E 	:	energy  [eV]
      @type   mparam 	:	minos config. parameters
      @param  mparam 	:	minos configuration information/parameters
    
      @rtype	        :	float
      @return	        :	dPhi/dE [cm^-2 GeV^-1 POT^-1]
    """
    return nb.NUMIflux(E/mparam.GeV,mparam.neutype)*nb.FNratio(E/mparam.GeV)
    #return nb.NUMIflux(E/mparam.GeV,mparam.neutype)*(mparam.NearDet.distance/mparam.FarDet.distance)**2
    #return nb.NUMIflux_binned(E/mparam.GeV,mparam.neutype)*(mparam.NearDet.distance/mparam.FarDet.distance)**2   
  
def NevtSpectraMINOS(E,mparam):
    """ Calculates the event number spectra at MINOS far detector.
  
      @type   E 	:	float
      @param  E 	:	energy  [eV]
      @type   mparam 	:	minos config. parameters
      @param  mparam 	:	minos configuration information/parameters
    
      @rtype	        :	float
      @return	        :	dN/dE   [event GeV^-1]
    """
    oscconf = mparam.NeuOscConf
    
    if mparam.neutype == 0:
        POT = mparam.Beam.POT_neutrino
    elif mparam.neutype == 1:
        POT = mparam.Beam.POT_antineutrino
    else :
        print "Wrong neutrino type."
        quit()
        
    if oscconf.oscillation :
        if oscconf.matter :
            body = bd.Earth()            
        else :
            body = bd.Vacuum()
            
        track = body.track(mparam.baseline*mparam.km,mparam.baseline*mparam.km,mparam.baseline*mparam.km)        
        if oscconf.generations :
            return FluxAtFarDetector(E,mparam)*(TargetNumber(mparam)*xs.nuMINOSxsection_CC(E/mparam.GeV,mparam.neutype))*POT*oa.NeuOsc3g(1,1,E,track,body,mparam)*MINOSEfficiency(E,mparam)
        else :
            return FluxAtFarDetector(E,mparam)*(TargetNumber(mparam)*xs.nuMINOSxsection_CC(E/mparam.GeV,mparam.neutype))*POT*oa.NeuOsc2g(1,1,E,track,body,mparam)*MINOSEfficiency(E,mparam)
    else : 
        return FluxAtFarDetector(E,mparam)*(TargetNumber(mparam)*xs.nuMINOSxsection_CC(E/mparam.GeV,mparam.neutype))*POT*MINOSEfficiency(E,mparam)

  
def NevtMINOS(Emin,Emax,mparam):
    """ Calculates the event number between Emin and Emax at MINOS far detector.
    
      @type   Emin 	:	float
      @param  Emin 	:	energy  [eV]
      @type   Emax 	:	float
      @param  Emax 	:	energy  [eV]
      @type   mparam 	:	minos config. parameters
      @param  mparam 	:	minos configuration information/parameters
    
      @rtype	        :	float
      @return	        :	N       [event number]
    """
    
    return (integrate.quad(lambda E : NevtSpectraMINOS(E,mparam)*mparam.GeV**-1,Emin,Emax, limit = 400)[0])
        
if __name__ == '__main__':
    ## TESTING AREA ##
    import neutrinocommon.physconst.physicsconstants as PC
    import neutrinocommon.exp.minos.minosclass as mc
    
    minos = mc.MinosExperiment()
    minos.Beam.POT_neutrino = 1.0e17
    pc = PC.PhysicsConstants()
    Emin = minos.minenergy*pc.GeV
    Emax = minos.maxenergy*pc.GeV
    
    
