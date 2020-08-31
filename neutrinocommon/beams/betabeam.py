""" 
Author  : M.Bustamante
Date    : 11/FEB/2012

Implements a beta beam. Multiple materials are available.

Log :
- Modified 12/FEB/2012 by C.Arguelles.
    + Python implementation of Mauricio's code.
- Modified 20/MAR/2012 by C.Arguelles
    + Changed FluxAtDetector to return the neutrino flux
    intead of the transversed neutrino flux. Reviewed 
    compatibility with Mauricio old code. 
    
    + From now on to check or review this code refer to
    
    Neutrino nucleus interaction rates at a low-energy beta beam facility
    Julien Serreau, Cristina Volpe
    Phys.Rev. C70 (2004) 055502
    arXiv : hep-ph/0403293
    
"""   

import numpy as np
import scipy as sp
import scipy.interpolate as interpolate
import neutrinocommon.physconst.physicsconstants as PC

pc = PC.PhysicsConstants()

def FermiFunction(W,Z=2):
    """ Returns the Fermi function F(+/-Z,W) for Z = 2, evaluated using the	#
    method in the reference, correct to within 1%.			#
    
    Ref.: J. Phys. G: Nucl. Phys. 11, 359 (1985)			  

    @type   W  :  float
    @param  W  :  energy [eV]

    @rtype         :   float
    @return        :   fermi function value [dimensionless?]
    """    
    me = pc.electron_mass*pc.GeV
    alp = pc.alpha

    p = np.sqrt(W*W-me*me)
    a = 2.0*pc.PI*alp*Z
    b = a / (1.0-np.exp(-a))
    C = b-a
    d = 0.5*(b-1.0)
    return a*W/p + C/(1.0+d/(p*p))
    
class BetaBeam():
    name = "Beta beam"
    
    def SetInyectionRate(self,material = 'He',inyection_rate = None):
        """ Sets and returns the inyection rate
        corresponding to the given material. If
        an inyection rate is specified this value 
        will be used.
        
        Inyections rate are referenced here :
        
        Ref : PRC 70, 055502 (2004)

        @type   material  :  string
        @param  material  :  He : Helium, Ne: Neon
        @type   inyection_rate      :   float
        @param  inyection_rate      :   inyection rate [s^-1]        

        @rtype         :   float
        @return        :   inyection rate [s^-1]
        """
        if inyection_rate != None : 
            self.inyection_rate = inyection_rate/pc.sec
        else :
            if material == 'He':
                inyection_rate = 2.0e13/pc.sec
            elif material == 'Ne':
                inyection_rate = 8.0e11/pc.sec
            else :
                print "NC:BEAMS:BETABEAM:ERROR: Inyection rate for material '"+material+"' not found."
                quit()
                
            self.inyection_rate = inyection_rate
        
        return inyection_rate*pc.sec
        
    def SetBetaBeamGeometry(self,total_length = 1885.0,straight_length = 678.0):
        """ Sets the beta beam geometry. The 
        beta beam source geometry is define by its total
        length and the straight sections. Here we set those
        two in meters.        

        @type   total_length  :  float
        @param  total_length  :  Total source length [m].
        @type   straight_length  :  float
        @param  straight_length  :  Straight sections [m].

        @rtype         :   none
        @return        :   nothing.
        """    
        
        self.total_length = total_length*pc.meter
        self.straight_length = straight_length*pc.meter
        
    def SetMaterialDecayProperties(self,material = 'He'):
        """ Sets the material decay properties.        

        @type   material  :  string
        @param  material  :  He : Helium, Ne: Neon

        @rtype         :   none
        @return        :   nothing.
        """    
        if material == 'He':
            # Data for production through He++(6,2) --> Li+++(6,3) + e- + nu_e_bar
            # Lifetime of parent nuclei He(6,2)
            tau = 1.1638*pc.sec # [s]
            # Q value of the reaction
            Q = 3.5078*pc.MeV # [MeV]
            # Comparative half-life or ft value
            ft = 806.7*pc.sec # [s]
        elif material == 'Ne':
            # Data for production through Ne++(18,10) --> F-(18,9) + e+ + nu_e
            # Lifetime of parent nuclei Ne(18,9)
            tau = 1.6728*pc.sec # [s]
            # Q value of the reaction
            Q = 4.4344*pc.MeV # [MeV]
            # Comparative half-life or ft value
            ft = 1672.8*pc.sec # [s]
        else :
            print "NC:BEAMS:BETABEAM:ERROR: Inyection rate for material '"+material+"' not found."
            quit()
            
        self.tau = tau
        self.Q = Q
        self.ft = ft

        
    def FluxAtDetector(self,Enu,L,theta = 0.0):
        """ ReturFluxAtDetectorns the total neutrino flux through the detector, assuming that
        the detector is located far from the source (i.e. L>>source_length)
        
        Ref: - Phys. Rev. C 70, 055502 (2004)				
        
        Note : In Mauricio's Original code the transverse flux is given. Which is related to this flux by S, where S is the detector cross section. 

        @type   Enu  :  float
        @param  Enu  :  neutrino energy [eV]
        @type   L    :  float
        @param  L    :  distance from the ring to the detector [eV^-1]
        @type   theta:  float
        @param  theta:  angle of emission with respect to beam axis [rad].

        @rtype       :   float
        @return      :   flux [eV cm^-2 s^-1 in nat. units]
        """
        PI = pc.PI
            
        lt = self.total_length
        ls = self.straight_length
        
        iny = self.inyection_rate
        tau = self.tau
        
        # tau is the ion lifetine. This means that 
        # iny*tau is the number of ions in the ring

        return iny*tau*self.FluxInLabFrame(Enu,theta = theta)*(ls/lt)*(1/(4.0*PI*L*L))
        
    def FluxInLabFrame(self,Enu,theta = 0.0):
        """ Returns the neutrino flux in the lab frame.
        
        Ref: - Phys. Rev. C 70, 055502 (2004)				

        @type   Enu  :  float
        @param  Enu  :  neutrino energy [eV]
        @type   theta:  float
        @param  theta:  angle of emission with respect to beam axis [rad].

        @rtype       :   float
        @return      :   flux [eV/s in nat. units]
        """     
        gamma = self.gamma

        beta = np.sqrt(1.0-1.0/(gamma*gamma))
        f = gamma*(1.0-beta*np.cos(theta))
        return self.FluxInCMFrame(Enu*f)/f
    
    def FluxInCMFrame(self,Enu):
        """ Returns the neutrino flux in the rest frame.
        
        Ref: - Phys. Rev. C 70, 055502 (2004)				
             - hep-ph/0702175	

        @type   Enu  :  float
        @param  Enu  :  neutrino energy [eV]

        @rtype       :   float
        @return      :   flux [eV/s in nat. units]
        """           
        Ee = self.Q - Enu
        em = pc.electron_mass*pc.GeV

        if Ee >= em:
            b = pc.ln2 / em**5 / self.ft
            return b * Enu*Enu * Ee * np.sqrt(Ee*Ee-em*em) * FermiFunction(Ee)
        else:
            return 0.0
    
    def __init__(self,material = 'He', gamma = 25.0):
        """ Initialize beta beam.

        @type   material    :   string
        @param  material    :  He : Helium, Ne: Neon
        @type   gamma       :  float
        @param  gamma       :  neutrino lorentz boost factor [dimensionless]

        @rtype       :   none
        @return      :   nothing
        """
        
        self.material = material
        self.gamma = gamma
        
        self.SetInyectionRate(material)
        self.SetMaterialDecayProperties(material)
        self.SetBetaBeamGeometry()
        
if __name__ == '__main__':
    # EXAMPLE OF A BETABEAM

    BB = BetaBeam('He')
    
    E = 50.0*pc.MeV
    L = 1000.0*pc.km
    R = 4.0*np.sqrt(5.0)*pc.meter
    S = pc.PI*R**2
    
    print BB.FluxAtDetector(E,L)*(pc.MeV*pc.sec*pc.cm**2)
    print BB.FluxInLabFrame(E)*(pc.MeV*pc.sec*pc.cm**2)
    
    E = 80.0*pc.MeV
    
    print BB.FluxAtDetector(E,L)*(pc.MeV*pc.sec*pc.cm**2)    
    print BB.FluxInLabFrame(E)*(pc.MeV*pc.sec*pc.cm**2)
        
        
        
    
    

