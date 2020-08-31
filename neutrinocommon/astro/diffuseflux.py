""" 
Author  : C.A.Arguelles
Date    : 11/MAY/2011

This package contains definition of the neutrino diffuse flux compatible
with the IceCube experiment setup implemented in neutrinocommon.exp.icecube.icecube.

Log :
- Modified on 30/MAY/2012 by C.A.Arguelles
    + Fix comments format and added description.
    + Fix units to make it compatible with exp.icecube program.
"""

# python modules
import scipy.interpolate as interpolate
import numpy as np
# my modules
import neutrinocommon.physconst.physicsconstants as PC

# starting physconstant global variable for units change
pc = PC.PhysicsConstants()

class HalzenDiffuseFlux():
    """
    Implementation of the neutrino diffuse flux proposed by F. Halzen et al.
    in arXiv : 0802.0887.
    """
    
    def __init__(self,alpha,open_angle = 5.73225313576):
        """ Initializes the flux and sets the specral index.
        
        @type   alpha   :   float
        @param  alpha   :   spectral index         [dimensionless]
        @type   open_angle   :   float
        @param  open_angle   :   openning angle    [sr]        
        
        @rtype          :   None
        @return         :   None
        """
        self.open_angle = open_angle
        self.alpha = alpha
        
        return None
   
    def SpectrumNormalization(self,alpha,return_interpolator = False):
        """ Sets the spectrum normalization in natural units. Originally
        the normalization was given in [TeV^-1 cm^-2 s^-1]
        
        Description of the interpolator :
            Recieves    alpha               :   spectral index
            Returns     flux normalization  :   normalization
            
        Ref. arXiv : 0802.0887
            
        @type   alpha   :   float
        @param  alpha   :   spectral index         [dimensionless]
        @type   return_interpolator   :   boolean
        @param  return_interpolator   :  True : returns interpolator, False : returns interpolated value
    
        @rtype          :   float
        @return         :   normalization [TeV^-1 cm^-2 s^-1] in natural units.
        
        """        
        
        # data points extracted from Ref. Fig. 1
        alp = [1.9963553389,2.0218020786,2.0424771212,2.0631417603,2.0869967784,
               2.1092531392,2.1315095001,2.1553541147,2.1776104755,2.2030364085,
               2.2237079832,2.2507187022,2.2761446351,2.2999788464,2.3269860976,
               2.3524085627,2.3794158139,2.4064195973,2.4413646499,2.4715380054,
               2.4985383210,2.5207634716,2.5493485732,2.5826880330,2.6128509850,
               2.6430104693,2.6763464613,2.7176098513,2.7541154154,2.7985483776,
               2.8302822445,2.8667704696,2.8985043365,2.9318195217,2.9667194929]
        
        dNdE = [13.0328227571,12.5776805252,12.2100656455,11.8949671772,
                11.4748358862,11.124726477,10.7746170678,10.4070021882,
                10.056892779,9.7067833698,9.3566739606,9.0065645514,
                8.6564551422,8.341356674,8.0087527352,7.6761487965,
                7.3435448578,7.0284463895,6.6258205689,6.3107221007,
                6.0131291028,5.8205689278,5.52297593,5.2253829322,
                4.9628008753,4.7177242888,4.4376367615,4.1400437637,
                3.8599562363,3.5623632385,3.3698030635,3.1772428884,
                2.9846827134,2.8096280088,2.6345733042]
        
        # creating interpolating function
        inter = interpolate.interp1d(alp,dNdE)
        # defining unit
        unit = pc.TeV**-1*pc.cm**2*pc.sec**-1
        
        # evaluating interpolator
        if return_interpolator :
            return inter
        else : 
            if (alpha <= alp[0]):
                return inter(alp[0])*1.0e-12*unit
            elif (alpha >= alp[-1]):
                return inter(alp[-1])*1.0e-12*unit
            else:
                return inter(alpha)*1.0e-12*unit
   
    def DiffuseFlux(self,Enu):
        """ Calculates the neutrino diffuse flux described in reference
        
        arXiv : 0802.0887
        
        @type   Enu   :   float
        @param  Enu   :   neutrino energy [eV]

        @rtype          :   float
        @return         :   diffuse flux [TeV^-1 cm^-2 s^-1 sr^-1] in natural units
        
        """           
        open_angle = self.open_angle
        alpha = self.alpha
        phi = (self.SpectrumNormalization(alpha)*(Enu/pc.GeV)**(-alpha))/open_angle
        return phi
    
    ################ ICECUBE DETECTOR COMPATIBILITY FUNCTIONS ##################
    
    def numu_flux(self,Enu):
        """ Returns the muon neutrino flux at the detector.
        
        @type   Enu   :   float
        @param  Enu   :   neutrino energy [eV]

        @rtype          :   float
        @return         :   diffuse flux [TeV^-1 cm^-2 s^-1 sr^-1] in natural units
        
        """
        return self.DiffuseFlux(Enu)
        
if __name__ == '__main__':
    #omg = 2.0*np.pi*(1-np.cos(icecube.thmax))
    #print omg
    omg = 1.0
    unit = pc.GeV**-1*pc.cm**2*pc.sec**-1
    HDF = HalzenDiffuseFlux(2.0)
    print HDF.SpectrumNormalization(2.0)/unit
    #print HDF.SpectrumNormalization(3.0)/unit
    
    E = 100.0*pc.GeV
    print HDF.numu_flux(E)/unit
    
    import neutrinocommon.exp.icecube.icecube as ice
    T = 1.0*pc.year
    print ice.NevtSpectra(HDF,0,T,pc)
    
    
    
        