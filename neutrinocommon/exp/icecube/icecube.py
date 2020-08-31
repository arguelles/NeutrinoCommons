"""
Author  : C.A. Arguelles
Date    : 10/MAY/2011

Basic functions to calculate event numbers at IceCube.

Log:

- Modified on 10/FEB/2012 by C.Arguelles
    + Change references to cross section to make it compatible with new cross 
    section file naming.
- Modified on 23/ABR/2012 by C.Arguelles
    + Fixed linking to data files.
    + Fixed cross section units according to the new output
    See neutrinocommon.neu.xsections.
- Modified on 30/MAY/2012 by C.Arguelles
    + Included areas optimized for dark matter searches
    + Fixed comments in IceCube class definition.
"""

import scipy.interpolate as interpolate
import numpy as np
#my modules
import neutrinocommon
import neutrinocommon.tools.generaltools as gt
import neutrinocommon.neu.xsections as xs

# global variables
filepath = neutrinocommon.exp.icecube.__path__
global_datapath = filepath[0]+"/data/"

class Icecube():
    """
    Icecube class containing basic IceCube definitions.
    """
    # Main IceCube Class
    # Icecube parameters
    Enumin      = 1.0e2     # [GeV]
    Enumax      = 1.0e8     # [GeV]
    Ebinsize    = 10.0      # [GeV]
    Emumin      = 1.0e2     # [GeV]
    Volume      = 1.0e6     # [m^3]
    Area        = 1.0e6     # [m^2] # hand estimation a= 2.5*np.pi*1.0e6
    Depth       = 2.0e3     # [m]
    thmax       = 1.483     # [radians]
    T           = 1.0       # [years]
    
    # Ice density    
    rho_ice     = 0.918     # [gr/cm^3]
    
    # Data arrays
    EMK = []        # Earthmatterkernel
    ESF = []        # Earthshadowfactor
    # Program flags
    Datgen = True
    # Miscelania
    globalfact = 1.0
    
    def Aeffmu(self,E):
        """ Returns IceCube effective muon detection area at level 2 cuts to E^-2
    
        @type   E       :   float
        @param  E       :   muon energy         [eV]
        
        @rtype          :   float
        @return         :   IceCube effective area [m^2]
        """
        
        E_GeV = E/pc.GeV
        
        if(1.0e1<=E_GeV<1.0e2):
            return 0.2692*1.0e6
        elif(1.0e2<=E_GeV<1.0e3):
            return 0.6385*1.0e6
        elif(1.0e3<=E_GeV<1.0e4):
            return 0.8462*1.0e6
        elif(1.0e4<=E_GeV<1.0e5):
            return 1.0231*1.0e6
        elif(1.0e5<=E_GeV<1.0e6):
            return 1.1615*1.0e6
        elif(1.0e6<=E_GeV<1.0e7):
            return 1.3385*1.0e6
        elif(1.0e7<=E_GeV<=1.0e8):
            return 1.50*1.0e6
        else :
            print "NC:EXP:ICECUBE:ERROR:Aeffmu: Muon effective area out of range."
            quit()
            
    def Aeffmu_trigger(self,E,use_binned = True):
        """ Returns IceCube effective muon detection area at trigger level.
        
        Fig. 5
        Sensitivity of the IceCube Detector to
        Astrophysical Sources of High Energy Muon Neutrinos
        arXiv : astro-ph/0305196
    
        @type   E       :   float
        @param  E       :   muon energy         [eV]
        @type   use_binned       :   boolean
        @param  use_binned       :   True : binned function, False : interpolated function
        
        @rtype          :   float
        @return         :   IceCube effective area [m^2] in natural units
        """
        units = pc.cm**2
        
        if use_binned :
            return self.Aeffmu_trigger_binned(E/pc.GeV)*units
        else :
            return self.Aeffmu_trigger_interpolated(E/pc.GeV)*units
        
            
    def Aeffmu_trigger_binned(self,E):
        """ Returns IceCube binned effective muon detection area at trigger level.
        
        Fig. 5
        Sensitivity of the IceCube Detector to
        Astrophysical Sources of High Energy Muon Neutrinos
        arXiv : astro-ph/0305196
    
        @type   E       :   float
        @param  E       :   muon energy         [GeV]
        
        @rtype          :   float
        @return         :   IceCube effective area [m^2]
        """
        if(0.0 <=E< 1.0e1):
            return 0.4923*1.0e6
        elif(1.0e1<=E<1.0e2):
            return 1.4945*1.0e6
        elif(1.0e2<=E<1.0e3):
            return 2.5319*1.0e6
        elif(1.0e3<=E<1.0e4):
            return 3.5165*1.0e6
        elif(1.0e4<=E<1.0e5):
            return 4.5010*1.0e6
        elif(1.0e5<=E<1.0e6):
            return 5.4857*1.0e6
        elif(1.0e6<=E<1.0e7):
            return 6.4879*1.0e6
        elif(1.0e7<=E<=1.0e8):
            return 7.4901*1.0e6
        else :
            print "NC:EXP:ICECUBE:ERROR:Aeffmu_trigger:Muon effective area out of range."
            quit()
            
    def Aeffmu_trigger_interpolated(self,E):
        """ Returns IceCube interpolated effective muon detection area at trigger level.
    
        @type   E       :   float
        @param  E       :   muon energy         [GeV]
        
        @rtype          :   float
        @return         :   IceCube effective area [m^2]
        """
        Enu     = [0.0,0.5e1,1.5e1,1.5e2,1.5e3,1.5e4,1.5e5,1.5e6,1.5e7,1.5e8]
        Aeff    = [0.4923*1.0e6,0.4923*1.0e6,1.4945*1.0e6,2.5319*1.0e6,3.5165*1.0e6,4.5010*1.0e6,5.4857*1.0e6,6.4879*1.0e6,7.4901*1.0e6,7.4901*1.0e6]
        inter   = interpolate.interp1d(Enu,Aeff)
        return inter(E)
        
    def NevtSpectra(self,flux,neutype = 'all',event_direction = 'UP'):
        """ Returns the IceCube event number considering the trigger level
        effective muon areas.
    
        @type   flux       :   flux object
        @param  flux       :   contains a flux object with a numu_flux rutine
        @type   neutype    :   string
        @param  neutype    :   'neutrino': neutrino flux, 'antineutrino':antineutrino flux, 'all' : neutrino + antineutrino
        @type   event_direction :   string
        @param  event_direction :   'UP' : upgoing neutrinos, 'DOWN' : downgoing neutrinos, 'TOTAL' : downgoing + upgoing
        
        @rtype          :   float
        @return         :   event number
        """
        
        avg_flux = lambda nuflux,anuflux,xs,axs,range,arange : (nuflux*xs*range + anuflux*axs*arange)/(xs*range+axs*arange)
        
        if event_direction == 'UP':
            if neutype == 'all':
                phi = Aeffmu_trigger(E)*avg_flux(flux.numu_flux(E,'neutrino'),flux.numu_flux(E,'antineutrino'),
                                                 xs.nuDISxsection_CC_Tbl(E,0),xs.nuDISxsection_CC_Tbl(E,1),
                                                 MuonRangeEnu(E,0,param),MuonRangeEnu(E,1,param))
            elif neutype == 'neutrino':
                phi = Aeffmu_trigger(E)*avg_flux(flux.numu_flux(E,'neutrino'),0.0,
                                                 xs.nuDISxsection_CC_Tbl(E,0),xs.nuDISxsection_CC_Tbl(E,1),
                                                 MuonRangeEnu(E,0,param),MuonRangeEnu(E,1,param))
            elif neutype == 'antineutrino':
                phi = Aeffmu_trigger(E)*avg_flux(0.0,flux.numu_flux(E,'antineutrino'),
                                                 xs.nuDISxsection_CC_Tbl(E,0),xs.nuDISxsection_CC_Tbl(E,1),
                                                 MuonRangeEnu(E,0,param),MuonRangeEnu(E,1,param))
        elif event_direction == 'DOWN':
            if neutype == 'all':
                phi = Aeffmu_trigger(E)*avg_flux(flux.numu_flux(E,'neutrino'),flux.numu_flux(E,'antineutrino'),
                                                 xs.nuDISxsection_CC_Tbl(E,0),xs.nuDISxsection_CC_Tbl(E,1),
                                                 MuonRangeEnu(E,0,param),MuonRangeEnu(E,1,param))
            elif neutype == 'neutrino':
                phi = Aeffmu_trigger(E)*avg_flux(flux.numu_flux(E,'neutrino'),0.0,
                                                 xs.nuDISxsection_CC_Tbl(E,0),xs.nuDISxsection_CC_Tbl(E,1),
                                                 MuonRangeEnu(E,0,param),MuonRangeEnu(E,1,param))
            elif neutype == 'antineutrino':
                phi = Aeffmu_trigger(E)*avg_flux(0.0,flux.numu_flux(E,'antineutrino'),
                                                 xs.nuDISxsection_CC_Tbl(E,0),xs.nuDISxsection_CC_Tbl(E,1),
                                                 MuonRangeEnu(E,0,param),MuonRangeEnu(E,1,param))
        elif event_direction == 'TOTAL':
            if neutype == 'all':
                flux.numu_flux(E,'neutrino')+flux.numu_flux(E,'antineutrino')
            elif neutype == 'neutrino':
                flux.numu_flux(E,'neutrino')
            elif neutype == 'antineutrino':
                flux.numu_flux(E,'antineutrino')
        else :
            print "NC:EXP:ICECUBE:ICECUBE:NevtSpectra: Invalid event_direction " + event_direction
            quit()
        
        return phi

IC = Icecube()

def FullyContainedNevtSpectra(flux,neutype,T,param,use_effective_area = False, flux_in_natural_units = True):
    """ Calculates event number spectra for UP going neutrinos.
    
    Spectra:
    
    dPhi/dE_nu = T * Phi_nu * cross_section * nucleon_num
    
    where
    
    nucleon_num = (Aeff_mu * R_muon)*density*Na

    @type   flux    :   flux object
    @param  flux    :   flux.neuflux(E,neutype) returns flux in [GeV^-1 m^-2 s^-1] or [m^-2 s^-1]
    @type   flavor  :   integer
    @param  flavor  :   neutrino flavor (0:e,1:mu,2:tau)
    @type   T       :   float
    @param  T       :   exposure time [yr]
    @type   param   :   physicsconstants
    @param  param   :   set of physical parameters to be used.
    @type   use_effective_area   :   boolean
    @param  use_effective_area   :   set of physical parameters to be used.
    @type   flux_in_natural_units   :   boolean
    @param  flux_in_natural_units   :   set of physical parameters to be used.    
    
    @rtype          :   float
    @return         :   dPhi/dE_nu at Icecube or Phi(E) if total flux is given.
    """
    Enu 	= np.arange(flux.minenergy/param.GeV,flux.maxenergy/param.GeV,0.1)
    T 		= (param.yearstosec)*Icecube.T
    
    # Here E is in GeV, some conv. are made with param.GeV.
    evt_spectra = map(lambda E : T*flux.numu_flux(E*param.GeV,neutype)*(xs.nuDISxsection_CC_Tbl(E,neutype)/param.meter**2)*NucleonNumber(MuonEnergy(E*param.GeV,neutype,param),param),Enu)  
    
    return evt_spectra

def MuonRangeEnu(E,neutype,param):
    """ Calculates muon range on ice for a given neutrino energy.

    @type   E       :   float
    @param  E       :   neutrino energy [eV]
    @type   neutype :   integer
    @param  neutype :   0 : neutrino - 1 : antineutrino    
    @type   param   :   physicsconstants
    @param  param   :   set of physical parameters to be used.
    
    @rtype       :   float
    @return      :   range [eV^-1]
    """
    return MuonRange(MuonEnergy(E,neutype,param),param)*param.meter

def MuonEnergy(E,neutype,param):
    """ Average energy of the resulting muon from QE, RES and DIS \nu-N interactions.    

    @type   E       :   float
    @param  E       :   neutrino energy [eV]
    @type   neutype :   integer
    @param  neutype :   0 : neutrino - 1 : antineutrino
    @type   param   :   physicsconstants
    @param  param   :   set of physical parameters to be used.
    
    @rtype          :   float
    @return         :   muon energy [eV]
    """
    return (1.0-xs.ymuCC_avg(E/param.GeV,neutype))*E    
    
def MuonRange(E,param):
    """ Calculates muon range on ice using PDG data using the muon energy.
    
    Ref. PDG data for Ice
    http://pdg.lbl.gov/2010/AtomicNuclearProperties/HTML_PAGES/325.html

    @type   E       :   float
    @param  E       :   muon energy [eV]
    @type   param   :   physicsconstants
    @param  param   :   set of physical parameters to be used.
    
    @rtype       :   float
    @return      :   range [m]
    """
    # load muon range pdg-data #
    filename = "muonloss_325.dat"
    file     = open(global_datapath + filename,'r')
    h,dat    = gt.hreadfilev2(file)

    rho_ice  = IC.rho_ice*param.gr/param.cm**3
    dat = np.array(dat)
    E_mu = dat[:,0]
    Range = dat[:,-3]*param.gr/param.cm**2
    
    inter_range = interpolate.interp1d(E_mu,Range)
    
    return (inter_range(E/param.MeV)/rho_ice)/param.meter
    
def MuonRange_IC(Emu,param):
    """ Calculates muon range using IC energy loss coefficients.

    @type   Emu     :   float
    @param  Emu     :   muon energy [eV]
    @type   param   :   physicsconstants
    @param  param   :   set of physical parameters to be used.
    
    @rtype       :   float
    @return      :   range [m]
    """
    
    a = 2.0e-3*param.GeV*param.cm**-1
    b = 3.9e-6*param.cm**-1
        
    Emumin = 1.0*param.GeV#IC.Emumin*param.GeV
    
    return (1.0/b)*np.log((a+b*Emu)/(a+b*Emumin))/param.meter
    
def MuonRange_SM(Emu,param):
    """ Calculates muon range using IC energy loss coefficients from
    Smirnov paper.
    
    arXiv : 1104.1390

    @type   Emu     :   float
    @param  Emu     :   muon energy [eV]
    @type   param   :   physicsconstants
    @param  param   :   set of physical parameters to be used.
    
    @rtype       :   float
    @return      :   range [m]
    """
    
    a = 0.24*param.GeV*param.meter**-1
    b = 3.3e-4*param.meter**-1
        
    Emumin = 1.0*param.GeV#IC.Emumin*param.GeV
    
    return (1.0/b)*np.log((a+b*Emu)/(a+b*Emumin))/param.meter       
    
def NucleonNumber(E,param):
    """ Calculates the number of nucleons in an Icecube-like detector
    using Icecube effective muon area and the muon range.

    @type   E       :   float
    @param  E       :   muon energy         [eV]
    @type   param   :   physicsconstants
    @param  param   :   set of physical parameters to be used.
    
    @rtype          :   float
    @return         :   nucleon number [dimensionless]
    """
    vol_eff     = IC.Aeffmu_trigger_inter(E/param.GeV)*MuonRange(E,param)*param.meter**3
    density     = IC.rho_ice*param.gr/param.cm**3
    mass_eff    = vol_eff*density
    
    nucleon_num = (mass_eff/param.gr)*param.Na
    
    return nucleon_num

################################################################################
########################## DARK MATTER OPTIMIZED AREAS #########################
################################################################################

def EffectiveArea(E,channel,use_binned = True):
    """ Returns IceCube effective area for DM search. Reference
    
    arXiv : 1111.2738
    M. Danninger's contribution, fig.2

    @type   E       :   float
    @param  E       :   neutrino energy         [GeV]
    @type   channel :   string
    @param  channel :   'WW' or 'bb'
    @type   use_binned :   boolean
    @param  use_binned :   True: will use bin function, False: will use interpolated function
    
    @rtype          :   float
    @return         :   IceCube effective area  [cm^2]
    """
    
    if use_binned : 
        return EffectiveAreaBinned(E,channel)
    else :
        return EffectiveAreaInterpolated(E,channel)
    
def EffectiveAreaBinned(E,channel):
    """ Returns IceCube effective area binned for DM search. Reference
    
    arXiv : 1111.2738
    M. Danninger's contribution, fig.2

    @type   E       :   float
    @param  E       :   neutrino energy         [GeV]
    @type   channel :   string
    @param  channel :   'WW' or 'bb'
    
    @rtype          :   float
    @return         :   IceCube effective area [cm^2]
    """    

    E_low = [ 0., 20.,         28.5272,     40.9201,    58.0387,    82.7841,    118.747   ]
    A_low = [ 0., 2.55067e-05, 9.74786e-05, 0.00033068, 0.00112178, 0.00282499, 0.0040391 ]

    E_high = [0.,  20.,   28.2074,     40.4614,    58.0387,    82.7841,    118.08,
              168.425,    240.234,     342.661,    488.757,    701.085	]
    A_high = [0., 1.1508e-06, 1.7316e-05, 0.000143587, 0.000995753, 0.00482968, 0.0149831,
              0.033493,   0.0572606,  0.0978944,   0.153053,   0.253984 ]

    if channel == 'bb':
      j = 0
      while j < len(A_low)-1 and E > E_low[j+1]:
        j = j + 1
      return A_low[j]
    elif channel == 'WW':
      j = 0
      while j < len(A_high)-1 and E > E_high[j+1]:
        j = j + 1
      return A_high[j]
    else:
      print "NC:EXP:ICECUBE:ERROR:EffectiveAreaBinned: Invalid channel: ", channel
      return 0.0
    
    
def EffectiveAreaInterpolated(E,channel):
    """ Returns IceCube effective area binned for DM search. Reference
    
    arXiv : 1111.2738
    M. Danninger's contribution, fig.2

    @type   E       :   float
    @param  E       :   muon energy         [GeV]
    @type   channel :   string
    @param  channel :   'WW' or 'bb'
    
    @rtype          :   float
    @return         :   IceCube effective area [cm^2]
    """
    
    E_low = [ 0., 20.,         28.5272,     40.9201,    58.0387,    82.7841,    118.747   ]
    A_low = [ 0., 2.55067e-05, 9.74786e-05, 0.00033068, 0.00112178, 0.00282499, 0.0040391 ]

    E_high = [0.,  20.,   28.2074,     40.4614,    58.0387,    82.7841,    118.08,
              168.425,    240.234,     342.661,    488.757,    701.085	]
    A_high = [0., 1.1508e-06, 1.7316e-05, 0.000143587, 0.000995753, 0.00482968, 0.0149831,
              0.033493,   0.0572606,  0.0978944,   0.153053,   0.253984 ]
    
    inter_low = interpolate.interp1d(E_low,A_low)
    inter_high = interpolate.interp1d(E_high,A_high)
    
    if channel == 'bb':
        if E <= E_low[0]:
            return A_low[0]
        elif E >= E_low[-1]:
            return A_low[-1]
        else:
            return inter_low(E)
    elif channel == 'WW':
        if E <= E_high[0]:
            return A_high[0]
        elif E >= E_high[-1]:
            return A_high[-1]
        else:
            return inter_high(E)
    
if __name__ == '__main__':
    # TEST AREA #
    import physicsconstants as PC
    import matplotlib.pyplot as plt
    import DM 
    pc = PC.PhysicsConstants()
    E  = 100.0*pc.MeV
    
    #Enu = np.arange(1.0,100.0,0.5)
    #nuc_num = map(lambda EE : NucleonNumber(EE*pc.GeV,pc),Enu)
    #
    #plt.plot(Enu,nuc_num)
    #plt.show()
    ch  = 'tautau'
    DMm = 100.0*pc.GeV
    
    DM_flux = DM.DMFluxAtDetector(ch,DMm,pc)

    
    
    
    
    
    
    
    
    
    
