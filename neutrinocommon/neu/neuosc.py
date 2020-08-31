""" 
Author  : C.A. Arguelles
Date    : 10/MAY/2011

This package contains definition of the Hamiltonian 
and lepton mixing matrixes. It also contains functions
to calculated and store the neutrino oscillation probabilities
as well as survival probabilities using the RK method.

Note : Density matrix formalism not implemented. Needed for
decoherence calculations, precise estimations of CC/NC
interactions. Currently decoherence formulae is implemented
on netruinocommon.neu.oscana and CC/NC interactions can be 
estimated using a MC.

Log :
- Modified on 9/FEB/2012 by C.Arguelles
    + Added error handling in order to run without 
    the cython libraries.
- Modified on 23/ABR/2012 by C.Arguelles
    + Fix the cross section usage in order to be consistant
    with the new cross section functions output. 

    Modified :    
        * CalNeuOscGSL
        * GenDataSurProb
    
    See neutrinocommon.neu.xsections.    
"""
# python standard modules
import numpy as np
import scipy as sp
from scipy import integrate,interpolate
from scipy.linalg.matfuncs import expm2
import datetime
import os
import fnmatch
try :
    import pygsl.odeiv as odeiv
    PYGSL_ERROR = False
except ImportError:
    print "NC:NEU:NEUOSC:ERROR: Loading GSL python interface : pygsl."
    PYGSL_ERROR = True
# my modules
import neutrinocommon.tools.generaltools as gt
import neutrinocommon.physconst.physicsconstants as PC
import neutrinocommon.astro.body as bd
import neutrinocommon.neu.xsections as xs
try :
    # disable the nusigma interface. need to integrate better. 
    # import neutrinocommon.shareobjects.optneuosc as ono
    OPTNEUOSC_ERROR = True #False
except ImportError:
    print "NC:NEU:NEUOSC:ERROR: Loading optimize neutrino oscillation functions : optneuosc."
    OPTNEUOSC_ERROR = True
try : 
    # disable the tauola interface. need to integrate better. 
    # import neutrinocommon.shareobjects.taudecay as td
    TAUOLA_ERROR = True #False
except ImportError:
    print "NC:NEU:NEUOSC:ERROR: Loading TAUOLA interface : taudecay."
    TAUOLA_ERROR = True
try : 
    import neutrinocommon.shareobjects.nudsde as oxs
    NUSIGMA_ERROR = False
except ImportError:
    print "NC:NEU:NEUOSC:ERROR: Loading NUSIGMA interface : nudsde."
    NUSIGMA_ERROR = True

#===============================================================================#
# Classes                                                                       #
#===============================================================================#

class ParamsContainer():
    """
    This is an auxiliary class that its used in #L{CalNeuOscGSL}.
    """
    # basic storage
    track   = None
    param   = None
    flavorM2 = None
    body    = None
    E       = None
    
    # used for optimization
    H0 = None
    S  = None
    rpos = None
    cH = None
    cS = None

class NeuFlavor():
    """
    This class defines a neutrino of a given flavor in 3 generations.
    
    NOTE : This class is usually no longer used. Now the easier 0:e,1:muon,2:tau convention is used.
    """
    def __init__(self,e,mu,tau):
        # set initial flavor state
        # e, mu, tau complex numbers
        self.e = e
        self.mu = mu
        self.tau = tau
        
    def mixmass(self,mixmatrix):
        nu1 = self.e*mixmatrix.UI[0,0]+self.mu*mixmatrix.UI[0,1]+self.tau*mixmatrix.UI[0,2]
        nu2 = self.e*mixmatrix.UI[1,0]+self.mu*mixmatrix.UI[1,1]+self.tau*mixmatrix.UI[1,2]
        nu3 = self.e*mixmatrix.UI[2,0]+self.mu*mixmatrix.UI[2,1]+self.tau*mixmatrix.UI[2,2]
        return neumass(nu1,nu2,nu3)
    
    def vec(self,param):
        # column vector format
        neu = np.zeros([param.numneu,1],complex)
        neu[0,0] = complex(self.e,0.0)
        neu[1,0] = complex(self.mu,0.0)
        neu[2,0] = complex(self.tau,0.0)    
        return neu
    
class NeuMass():
    """
    This class defines a neutrino of a given mass state in 3 generations.
    """
    def __init__(self,nu1,nu2,nu3):
        # set initial flavor state
        # e, mu, tau complex numbers
        self.nu1 = nu1
        self.nu2 = nu2
        self.nu3 = nu3
        
    def mixflavor(self,mixmatrix):
        e = self.nu1*mixmatrix.U[0,0]+self.nu2*mixmatrix.U[0,1]+self.nu3*mixmatrix.U[0,2]
        mu = self.nu1*mixmatrix.U[1,0]+self.nu2*mixmatrix.U[1,1]+self.nu3*mixmatrix.U[1,2]
        tau = self.nu1*mixmatrix.U[2,0]+self.nu2*mixmatrix.U[2,1]+self.nu3*mixmatrix.U[2,2]
        return neuflavor(e,mu,tau)
    
    def vec(self,param):
        # column vector format
        neu = np.zeros([param.numneu,1],complex)
        neu[0,0] = complex(self.nu1,0.0)
        neu[1,0] = complex(self.nu2,0.0)
        neu[2,0] = complex(self.nu3,0.0)    
        return neu
        
class MixMatrix():
    """
    This class contains all the definition of the rotation matrix and lepton mixing parametrization
    to generate the lepton mixing matrix. An arbitrary number of neutrino flavors can be implemented.
    """
    def __init__(self,param):
        """ Initializing class.
        Basic initialization.
        @type   param   :   PhysicsConstants
        @param  param   :   set of physical parameters to be used.
    
        @rtype          :   MixMatrix
        @return         :   returns mixing matrix class.  
        """        
        if param.U == None:
            # if no mixing matrix is specified calculate it using the mixing angles.
            self.calcU(param)
            self.calcUCT(param)
        else:
            # checking mixing matrix dimension
            U_array = np.array(param.U)
            shape = U_array.shape
            x_dim = shape[0]
            y_dim = shape[1]            
            
            if x_dim == y_dim and x_dim == param.numneu :
                # dim are ok
                print "NC:NEU:NEUOSC:MixMatrix:WARNING:Using specified mixing matrix."
                self.U = U_array
                self.UCT = self.U.conjugate().T
            else : 
                print "NC:NEU:NEUOSC:MixMatrix:ERROR:Input mixing matrix dimensions and neutrino numbers are not compatible."
                quit()
      
    def R(self,i,j,cp,param):
        """ Rotation Matrix
        Calculates the R_ij rotations. Also incorporates CP-phases when necesary.
        @type   i       :   int
        @param  i       :   i-column.
        @type   j       :   int
        @param  j       :   j-row.
        @type   cp      :   int
        @param  cp      :   if cp = 0 : no CP-phase. else CP-phase = CP_array[cp]
    
        @rtype          :   numpy.array
        @return         :   returns the R_ij rotation matrix.
        """
        # if cp = 0 -> no complex phase
        # R_ij, i<j
        if(j<i):
            # no funny business
            l = i
            i = j
            j = l
        
        # rotation matrix generator
        R = np.zeros([param.numneu,param.numneu],complex)
        # diagonal terms
        for k in np.arange(0,param.numneu,1):
            if(k != i-1 and k != j-1):
                R[k,k] = 1.0
            else :
                R[k,k] = param.c[i,j]
        # non-diagonal terms
        if(cp != 0):
            sd = np.sin(param.dcp[cp])
            cd = np.cos(param.dcp[cp])
            faseCP = complex(cd,sd)
        else:
            faseCP = complex(1.0,0.0)
        
        R[i-1,j-1] = param.s[i,j]*faseCP.conjugate()
        R[j-1,i-1] = -param.s[i,j]*faseCP
        return R
        
    #==========================================================================#
    # ###################### BEGIN MIX MATRIX CALCULATION #################### #   
    #==========================================================================#
        
    def calcU(self,param):
        """ Defining the mixing matrix parametrization.
        @type   param   :   PhysicsConstants
        @param  param   :   set of physical parameters to be used.
    
        @rtype          :   None
        @return         :   Sets mixing matrix.
        """
        if(param.numneu == 3):
            self.U =  np.dot(self.R(2,3,0,param),np.dot(self.R(1,3,1,param),self.R(1,2,0,param)))
        elif(param.numneu == 4):
            self.U =  np.dot(self.R(3,4,0,param),np.dot(self.R(2,4,2,param),np.dot(self.R(1,4,0,param),np.dot(self.R(2,3,0,param),np.dot(self.R(1,3,1,param),self.R(1,2,0,param))))))
        elif(param.numneu == 5):
            self.U =  np.dot(self.R(4,5,0,param),np.dot(self.R(3,5,0,param),np.dot(self.R(2,5,0,param),np.dot(self.R(1,5,3,param),np.dot(self.R(3,4,0,param),np.dot(self.R(2,4,0,param),np.dot(self.R(1,4,2,param),np.dot(self.R(2,3,0,param),np.dot(self.R(1,3,1,param),self.R(1,2,0,param))))))))))
        elif(param.numneu == 6):
            # 3+3 twin-sterile-neutrino model
            self.U =  np.dot(self.R(3,6,0,param),np.dot(self.R(2,5,0,param),np.dot(self.R(1,4,0,param),np.dot(self.R(2,3,0,param),np.dot(self.R(1,3,1,param),self.R(1,2,0,param))))))
        else:
            print "Sorry, too many neutrinos. Not yet implemented! =(."
            quit()
            
        # antineutrino case
        if param.neutype == "antineutrino" :
            self.U = self.U.conjugate()
    
    #===========================================================================#
    # ######################## END MIX MATRIX CALCULATION ##################### #
    #===========================================================================#
        
    def calcUI(self,param):
        """ Inverse mixing matrix.
        @type   param   :   PhysicsConstants
        @param  param   :   set of physical parameters to be used.
    
        @rtype          :   None
        @return         :   Sets inverse mixing matrix.
        """    
        UU = np.matrix(self.U)
        UUI = UU.I
        self.UI = UUI
        
    def calcUT(self,param):
        """ Transpose mixing matrix.
        @type   param   :   PhysicsConstants
        @param  param   :   set of physical parameters to be used.
    
        @rtype          :   None
        @return         :   Sets transpose mixing matrix.
        """        
        self.UT = self.U.T
        
    def calcUCT(self,param):
        """ Dagger mixing matrix (=transpose + conjugate).
        @type   param   :   PhysicsConstants
        @param  param   :   set of physical parameters to be used.
    
        @rtype          :   None
        @return         :   Sets dagger mixing matrix.
        """    
        self.UCT = self.U.conjugate().T
        
    def calcRho(self,param):
        """ Calculates the density matrix
        
        rho = U \dagger U
	
	WARNING:
	This function doesn't make sense since it always
	return the unit matrix. Something wrong here.
        
        @type   param   :   PhysicsConstants
        @param  param   :   set of physical parameters to be used.
    
        @rtype          :   None
        @return         :   Sets dagger mixing matrix.
        """        
        self.Rho = np.dot(self.U,self.UCT)
        
#===============================================================================#
# Calculating the Osc. Probaility                                               #
#===============================================================================#        
        
def rvec(vv,param):
    """ Converts a complex n-vector into a real 2n-vector.
    
    Note : param.numneu sets the vector size.
    
    @type   vv      :   List
    @param  vv      :   complex vector.
    @type   param   :   PhysicsConstants
    @param  param   :   set of physical parameters to be used.

    @rtype          :   List
    @return         :   real vector.
    """ 
    neu = np.zeros([2*param.numneu],float)
    for i in np.arange(0,2*param.numneu,2):
        neu[i] = vv[i/2,0].real
        neu[i+1] = vv[i/2,0].imag
    return neu
    
def cvec(vv,param):
    """ Converts a real 2n-vector into a complex n-vector.
    
    Note : param.numneu sets the vector size.
    
    @type   vv      :   List
    @param  vv      :   complex vector.
    @type   param   :   PhysicsConstants
    @param  param   :   set of physical parameters to be used.

    @rtype          :   List
    @return         :   real vector.
    """ 
    neu = np.zeros([param.numneu,1],complex)
    for i in np.arange(0,2*param.numneu,2):
        neu[i/2,0] = complex(vv[i],vv[i+1])
    return neu           
       
def massM2(param):
    """ Mass term in the neutrino mass basis.
    
    @type   param   :   PhysicsConstants
    @param  param   :   set of physical parameters to be used.

    @rtype          :   numpy array
    @return         :   mass matrix in mass basis.
    """ 
    M2 = np.zeros([param.numneu,param.numneu],complex)
    for k in np.arange(1,param.numneu,1):
        M2[k,k] = param.dmsq[k+1]
    return M2

def flavorM2(param):
    """ Mass term in the neutrino flavor basis.
    
    @type   param   :   PhysicsConstants
    @param  param   :   set of physical parameters to be used.

    @rtype          :   numpy array
    @return         :   mass matrix in flavor basis.
    """ 
    mixing_matrix = MixMatrix(param)
    return np.dot(mixing_matrix.U,np.dot(massM2(param),mixing_matrix.UCT))
    
def flavorAcc(param,E,body,track):
    """ Neutrino matter potential.
    
    Note : Neutrino energy is not used in this function since the neutrino energy has
    been factorized an is on the hamiltonian.
    
    @type   param   :   PhysicsConstants
    @param  param   :   set of physical parameters to be used.
    @type   E       :   Float
    @param  E       :   neutrino energy [eV].
    @type   body    :   Body
    @param  body    :   body where the neutrino will propagate (e.g. Earth,Sun,...)
    @type   track   :   Track
    @param  track   :   trayectory within the body.

    @rtype          :   numpy array
    @return         :   neutrino matter potential matrix.
    """ 
   
    Acc = np.zeros([param.numneu,param.numneu],complex)
    ye = body.ye(track)
    
    CC = param.sqr2*param.GF*param.Na*param.cm**-3*body.density(track)*ye
    NC = CC*(-0.5*(1.0-ye)/ye)
    
    Acc[0,0] = complex(CC+NC,0.0)
    Acc[1,1] = complex(NC,0.0)
    Acc[2,2] = complex(NC,0.0)        
    
    if param.neutype == "neutrino":
        return Acc
    elif param.neutype == "antineutrino":
        return -Acc
    else :
        print "Wrong neutrino type."
        quit()
    
def flavorH(track,param,flavorM2,body,E):
    """ Hamiltonian in the flavor basis multiplied by i.

    Note : This is really the RHS of the Mikheev-Smirnov-Wolfenstein equation.

    @type   track   :   Track
    @param  track   :   trayectory within the body.    
    @type   param   :   PhysicsConstants
    @param  param   :   set of physical parameters to be used.
    @type   flavorM2:   numpy array
    @param  flavorM2:   Mass square difference in the flavor basis.
    @type   body    :   Body
    @param  body    :   body where the neutrino will propagate (e.g. Earth,Sun,...)
    @type   E       :   Float
    @param  E       :   neutrino energy [eV].    

    @rtype          :   numpy array
    @return         :   neutrino matter potential matrix.
    """ 
    ii = complex(0.0,-1.0)
    return (ii)*(flavorM2/(2.0*E)+flavorAcc(param,E,body,track))

def RHS(neuvec,x,track,param,flavorM2,body,E):
    """ Evaluates the RHS of Schrodinger's equation.

    @type   neuvec  :   List
    @param  neuvec  :   real neutrino state 2n-vector.
    @type   x       :   Float
    @param  x       :   position in the trayectory [eV^-1]
    @type   track   :   Track
    @param  track   :   trayectory within the body.        
    @type   param   :   PhysicsConstants
    @param  param   :   set of physical parameters to be used.
    @type   flavorM2:   numpy array
    @param  flavorM2:   Mass square difference in the flavor basis.
    @type   body    :   Body
    @param  body    :   body where the neutrino will propagate (e.g. Earth,Sun,...)
    @type   E       :   Float
    @param  E       :   neutrino energy [eV].    

    @rtype          :   List
    @return         :   real neutrino state 2n-vector
    """ 
    # Update the position in the trayectory
    track.x = x
    # Convert to the complex plane. And then evaluate the RHS.
    if not OPTNEUOSC_ERROR: 
        neuvec = ono.cvec(neuvec,param)
    else :  
        neuvec = cvec(neuvec,param) 
    rhs = np.dot(flavorH(track,param,flavorM2,body,E),neuvec)
    # Return to the real plane
    tmprhs = (np.array([[val[0].real,val[0].imag] for val in rhs])).flatten()       
    
    return tmprhs
    
def RHS_INT_GSL(x,neuvec,params):
    """ Evaluates the RHS of Schrodinger's equation.

    @type   x       :   Float
    @param  x       :   position in the trayectory [eV^-1]
    @type   neuvec  :   List
    @param  neuvec  :   real neutrino state 2n-vector.    
    @type   params  :   ParamsContainer
    @param  params  :   ParamsContainer
    
    @rtype          :   List
    @return         :   real neutrino state 2n-vector
    """ 
    ii = complex(0.0,1.0)
    track,param,flavorM2,body,E,H0 = params.track,params.param,params.flavorM2,params.body,params.E,params.H0
    
    # Right hand side in Schrodinger's equations in the interaction basis
    
    # moving around 
    x_prev      = params.rpos
    x_current   = x
    track.x = x_current
    
    ## this should be done in complex plane -- TMP
    neuvec = cvec(neuvec,param)
    
    # Decompose the Hamiltonian H = H0 + H1
    # where :
    #       H0 = H_(r_i)
    #       H1 = H_(r_i+1) - H_(r_i)
    
    H_current = ii*flavorH(track,param,flavorM2,body,E)
    H1 = H_current - H0

    # opt -> store S
    #S = S#expm2(-ii*H0*(x_current - x_prev))
    S = expm2(-ii*H0*(x_current - x_prev))

    interaction_H = np.dot(S,np.dot(H1,S.conjugate().T))
    
    rhs_int = (-ii)*interaction_H
    
    rhs = np.dot(rhs_int,neuvec)    
    tmprhs = (np.array([[val[0].real,val[0].imag] for val in rhs])).flatten()
    
    # storing stuff
    params.cH = H_current
    params.cS = S
    return tmprhs    
    
def OscProb(cvec,param):
    """ Calculates probability from a complex state vector.
    
    Note : param.numneu sets the vector size.    
    
    @type   cvec    :   List
    @param  cvec    :   complex vector.
    @type   param   :   PhysicsConstants
    @param  param   :   set of physical parameters to be used.

    @rtype          :   List
    @return         :   real probability vector.
    """ 
    p = [ cvec[i].real**2+cvec[i].imag**2 for i in range(param.numneu)]
    return p
    
def CompleteDecoherenceOscProb(fneu_flavor,PMNS,param):
    """ Calculates probability from a complex state vector on the 
    complete decoherence case.
    
    Note : param.numneu sets the vector size.    
    
    @type   fneu_flavor    :   List
    @param  fneu_flavor    :   complex vector.
    @type   PMNS    :   MixingMatrix
    @param  PMNS    :   lepton mixing matrix.
    @type   param   :   PhysicsConstants
    @param  param   :   set of physical parameters to be used.

    @rtype          :   List
    @return         :   real probability vector.
    """
    fneu_mass = np.dot(PMNS.UCT,fneu_flavor)
    prob_flv = []
    for flv in range(param.numneu):
        vfinal = PMNS.U[flv,:]
        prob = 0.0
        for i in range(param.numneu):
            prob = prob + gt.norm(fneu_mass[i])**2*gt.norm(vfinal[i])**2
        prob_flv.append(prob[0])
    return prob_flv    
    
def EvolveNeuState(ineu,E,track,body,flavorM2,param):
    """ Evolves the neutrino state following the trayectory indicated in track.
    The initial neutrino state is specified by ineu and its energy is E.

    @type   ineu    :   int
    @param  ineu    :   0:electron,1:muon,2:tau
    @type   E       :   float
    @param  E       :   neutrino energy [eV].       
    @type   track   :   Track
    @param  track   :   trayectory within the body.       
    @type   body    :   Body
    @param  body    :   body where the neutrino will propagate (e.g. Earth,Sun,...)    
    @type   flavorM2:   numpy array
    @param  flavorM2:   Mass square difference in the flavor basis.
    @type   param   :   PhysicsConstants
    @param  param   :   set of physical parameters to be used.

    @rtype          :   List
    @return         :   complex neutrino state n-vector
    """ 
    if (ineu == 0 ):
        ineu = NeuFlavor(1.0,0.0,0.0)
    elif (ineu == 1 ):
        ineu = NeuFlavor(0.0,1.0,0.0)
    elif (ineu == 2 ):
        ineu = NeuFlavor(0.0,0.0,1.0)
    
    if not OPTNEUOSC_ERROR : 
        fneu = cvec(integrate.odeint(ono.RHS,rvec(ineu.vec(param),param),[track.xini,track.xend],args = (track,param,flavorM2,body,E),mxstep = 30000000,atol = 1.0e-8)[1],param)
    else :
        fneu = cvec(integrate.odeint(RHS,rvec(ineu.vec(param),param),[track.xini,track.xend],args = (track,param,flavorM2,body,E),mxstep = 30000000,atol = 1.0e-8)[1],param) 
    return fneu    

def CalNeuOscPy(ineu,E,track,body,flavorM2,param,complete_decoherence = False):
    """ Calculates the probability vector [P_ie,P_imu,P_itau,P_is,...] where i is 
    the neutrino flavor specified by ineu.
    
    @type   ineu    :   int
    @param  ineu    :   flavor of the initial neutrino (0 : e, 1 : mu, 2 : tau)
    @type   E       :   float
    @param  E       :   neutrino energy [eV]
    @type   track   :   Track
    @param  track   :   trayectory within the body.       
    @type   body    :   Body
    @param  body    :   body where the neutrino will propagate (e.g. Earth,Sun,...)    
    @type   flavorM2:   numpy array
    @param  flavorM2:   Mass square difference in the flavor basis.
    @type   param   :   PhysicsConstants
    @param  param   :   set of physical parameters to be used.    
    @type   complete_decoherence   :   Boolean
    @param  complete_decoherence   :   If True will calculate probabilities assuming complete decoherence.
    
    @rtype          :   List
    @return         :   oscillation probability list.
    """
    # evolves the initial neutrino state
    frkneu = EvolveNeuState(ineu,E,track,body,flavorM2,param)
    # load the mixing matrix
    PMNS = MixMatrix(param)
    prob_flv = []
    
    if complete_decoherence :
        prob_flv = CompleteDecoherenceOscProb(frkneu,PMNS,param)
    else :
        #prob_flv = [ gt.norm(psi)**2 for psi in frkneu ]
        prob_flv = OscProb(frkneu,param)
        
    return map(float,prob_flv)
    
def CalNeuOscGSL(ineu,E,track,body,flavorM2,param,file = None,family = 0,NCCC_int = False,complete_decoherence = False,survival_prob = False,E_threshold = 1.0,optimization = True,only_osc_prob = True,abs_error = 1.0e-7, rel_error = 1.0e-6):
    """ Calculates the probability vector [P_ie,P_imu,P_itau,P_is,...] where i is 
    the neutrino flavor specified by ineu. Also returns the survival probability.
    
    Calculates the probabilities using the RK algorithm from the GSL library. When
    the optimization flag is True it uses rediagonalization at each step in order
    to increase speed. This technique is documented on C.A.Arguelles,J.Kopp (arXiv:1202.3431).
    
    When the NCCC_int flag is True neutrino interactions with matter will be considered,
    in such a way that if the neutrino experiences a NC-interaction its energy will
    be shifted accordingly, and if the neutrino experiences a CC-interaction in the
    electron or muon flavor state the event will be discarted. On the other hand, if
    the neutrino experiences a CC-interaction in the tau flavor tau regeneration is 
    implemented by calling this routine recursively.
    
    NC and CC interactions are implemented using the differential cross section 
    provided by NUSIGMA, here implemented by the nudsde nusigma-interface. Futhermore
    tau regeneration is implemented using TAUOLA via the taudecay tauola-interface.
    
    When called by the #L{GenerateMCEvents} this provides an easy way of generating
    montecarlo events on the NCCC_int = True mode.
    
    UPDATED: 
    - Updated to also return the survival probabilities at each step when survival_prob = True.

    @type   ineu    :   int
    @param  ineu    :   initial neutrino state
    @type   E       :   float
    @param  E       :   neutrino energy [eV]
    @type   track   :   track object
    @param  track   :   asociated trayectory in body
    @type   body    :   body
    @param  body    :   body with the asociated density profile.
    @type   flavorM2:   array
    @param  flavorM2:   M^2 in flavor basis
    @type   param   :   physicsconstants
    @param  param   :   set of physical parameters to be used.
    @type   file    :   file object
    @param  file    :   open file where data will be printed, if None will return result.
    @type   family  :   integer
    @param  family  :   returns the same integer for all neutrinos in the same event.
    @type   NCCC_int:   boolean
    @param  NCCC_int:   toggle NC-CC interactions ON/OFF
    @type   complete_decoherence:   boolean
    @param  complete_decoherence:   toggle complete decoherence interactions ON/OFF    
    @type   survival_prob:   boolean
    @param  survival_prob:   save intermediate survival probabilities.
    @type   E_threshold:   float
    @param  E_threshold:   Lower neutrino energy threshold. [GeV]
    @type   optimization:   boolean
    @param  optimization:   toogle rediagonalization optimization ON/OFF.
    @type   only_osc_prob:   boolean
    @param  only_osc_prob:   will only print the oscillation probabilities.
    @type   abs_error:   float
    @param  abs_error:   set absolute error in the RK.
    @type   rel_error:   float
    @param  rel_error:   set relative error in the RK.
    
    @rtype          :   List
    @return         :   oscillation probability list.    
    """
    # define complex i
    ii = complex(0.0,1.0)

    # saving neutrino type
    if param.neutype == "neutrino":
        neu = 0
    elif param.neutype == "antineutrino":
        neu = 1
        
    # setup params aux. container
    params = ParamsContainer()
    params.track, params.param,params.flavorM2,params.body,params.E = track,param,flavorM2,body,E
    
    if not OPTNEUOSC_ERROR:
        params.H0 = ii*ono.flavorH(track,param,flavorM2,body,E)
    else : 
        params.H0 = ii*flavorH(track,param,flavorM2,body,E)
    params.S  = np.identity(param.numneu, complex)
    
    # setup GSL solver
    dim     = int(2*param.numneu)
    stepper = odeiv.step_rkf45
    if not OPTNEUOSC_ERROR:
        step    = stepper(dim,ono.RHS_INT_GSL,None,args = params)
    else :        
        step    = stepper(dim,RHS_INT_GSL,None,args = params)
    control = odeiv.control_y_new(step,abs_error,rel_error)
    evolve  = odeiv.evolve(step,control,dim)
    
    # initial parameters
    neu_state           = [0.0]*dim
    neu_state[ineu*2]   = 1.0             # initial state
    h = 10.0*param.km                     # initial step size
    r = track.xini
    # survial array
    survival_prob_list = []
    
    E_in    = params.E/param.GeV
    
    # survival probability
    p_survival = 1.0
    
    # starting interaction counters
    m = 0
    n = 0
    k = 0
    
    while r < track.xend :
        # advancing one RK step
        r_old = r
        params.rpos = r
        r,h,neu_state = evolve.apply(r,track.xend,h,neu_state)
        r_new = r
        
        r_adv = r_new-r_old
        
        # ZAZ beginning rediagonalization optimization ZAZ
        # CALCULATE exp(-i H0 (r-r0))
        # where H0 = H(r_current) = cH
        S = sp.linalg.matfuncs.expm2(-ii*params.cH*(r_adv))
        
        # saving this for next calculation
        params.H0 = params.cH
        # UNDO LAST TRANSFORMATION : psi = S psi
        if not OPTNEUOSC_ERROR:
            neu_state = ono.rvec(np.dot(S.conjugate().T,ono.cvec(neu_state,param)),param)
        else : 
            neu_state = rvec(np.dot(S.conjugate().T,cvec(neu_state,param)),param)

        if NCCC_int : 
            # CC-NC MONTECARLO
            E_GeV = params.E/param.GeV
            # pick my random point (E,p/E)
            E_MC = np.random.uniform(0.0,E_GeV,1)[0]
            p_MC = np.random.uniform(0.0,1.0/E_GeV,1)[0]
            
            # THIS LINES TURN OUT TO BE THE SAME AS Na = 1/AMU BY DEFINITION
            #pfac = np.abs(body.density(track)*param.Na*(r_adv/param.cm))
            pfac = np.abs((body.density(track)/param.atomic_mass_unit)*(r_adv/param.cm))
            
            # active flavors probabilities
            pe = neu_state[0]*neu_state[0]+neu_state[1]*neu_state[1]
            pm = neu_state[2]*neu_state[2]+neu_state[3]*neu_state[3]
            pt = neu_state[4]*neu_state[4]+neu_state[5]*neu_state[5]
            
            # probability that is an active neutrino
            p_active = pe+pm+pt
            
            # Figuring out if the neutrino will undergo an interaction 
            # by stacking the diferential cross sections and selecting 
            # a random point in the plane dSigma/dE versus E plane.
            
            # Note : converting from natural units to [cm^2 GeV^-1].
            
            dsde_NC     = xs.dsde_NC_opt(E_GeV,E_MC,param)/(param.cm**2/param.GeV)
            dsde_CC_e   = xs.dsde_CC_opt(E_GeV,E_MC,param,flv = 0)/(param.cm**2/param.GeV)
            # Assuming that the muon diff. cross section is similar to the electron
            dsde_CC_m   = dsde_CC_e
            # else we should use this line :
            #dsde_CC_mu = xs.dsde_CC_opt(E_GeV,E_MC,param,flv = 1)/(param.cm**2/param.GeV)
            dsde_CC_t   = xs.dsde_CC_opt(E_GeV,E_MC,param,flv = 2)/(param.cm**2/param.GeV)
            
            # new xsections
            # Note : converting from natural units to [cm^2].
                        
            sig_NC = xs.nuDISxsection_NC_NusigmaInt(E_GeV,param)/(param.cm**2)
            sig_CC_e = xs.nuDISxsection_CC_NusigmaInt(E_GeV,0,param)/(param.cm**2)
            sig_CC_m = xs.nuDISxsection_CC_NusigmaInt(E_GeV,1,param)/(param.cm**2)
            sig_CC_t = xs.nuDISxsection_CC_NusigmaInt(E_GeV,2,param)/(param.cm**2)
            
            sig_total_active = sig_NC*p_active + sig_CC_e*pe + sig_CC_m*pm +sig_CC_t*pt
            
            if p_MC <= (dsde_NC*p_active/sig_total_active)*(1.0-np.exp(-pfac*sig_total_active)):                 
                # changing energy due to NC interaction
                params.E = E_MC*param.GeV
                m = m + 1
            elif p_MC <= ((dsde_NC*p_active + dsde_CC_e*pe + dsde_CC_m*pm)/sig_total_active)*(1.0-np.exp(-pfac*sig_total_active)):                
                # discart event because of CC interation with e or muon.
                k = k + 1
                neu_state   = [0.0]*(dim)
                result      = [0.0]*(dim/2)
                break
            elif p_MC <= ((dsde_NC*p_active + dsde_CC_e*pe + dsde_CC_m*pm+dsde_CC_t*pt)/sig_total_active)*(1.0-np.exp(-pfac*sig_total_active)):                
                # decay tau
                if not PC.tauola_init:
                    # decay and initialize tauola
                    neutrino_energies = td.TaudecayNeutrinos(E_MC)
                else:
                    # decay 
                    neutrino_energies = td.TaudecayNeutrinos(E_MC,1)
                    
                PC.tauola_init = True
                if neutrino_energies != []:
                    for i in range(0,len(neutrino_energies),2):
                        params.track.xini = r
                        if int(neutrino_energies[i]>0):
                            param.neutype = "neutrino"
                            if int(neutrino_energies[i]) == 12:
                                # create electron-neutrino
                                CalNeuOscGSL(0,neutrino_energies[i+1]*param.GeV,params.track,params.body,flavorM2,param,file,family)
                            elif int(neutrino_energies[i]) == 14:
                                # create muon-neutrino
                                CalNeuOscGSL(1,neutrino_energies[i+1]*param.GeV,params.track,params.body,flavorM2,param,file,family)
                            elif int(neutrino_energies[i]) == 16:
                                # create tau-neutrino
                                CalNeuOscGSL(2,neutrino_energies[i+1]*param.GeV,params.track,params.body,flavorM2,param,file,family)
                        else:
                            param.neutype = "antineutrino"
                            if int(neutrino_energies[i]) == -12:
                                # create electron-antineutrino
                                CalNeuOscGSL(0,neutrino_energies[i+1]*param.GeV,params.track,params.body,flavorM2,param,file,family)
                            elif int(neutrino_energies[i]) == -14:
                                # create muon-antineutrino
                                CalNeuOscGSL(1,neutrino_energies[i+1]*param.GeV,params.track,params.body,flavorM2,param,file,family)
                            elif int(neutrino_energies[i]) == -16:
                                # create tau-antineutrino
                                CalNeuOscGSL(2,neutrino_energies[i+1]*param.GeV,params.track,params.body,flavorM2,param,file,family)
                # count tau
                n = n + 1
                break
            if E_GeV <= E_threshold :
                # discart event because of energy treshold
                result = [0.0]*(dim/2)
                break
        elif survival_prob:
            E_GeV = params.E/param.GeV
            # calculate probability survival
            pfac = np.abs((body.density(track)/param.atomic_mass_unit)*(r_adv/param.cm))
            
            pe = neu_state[0]*neu_state[0]+neu_state[1]*neu_state[1]
            pm = neu_state[2]*neu_state[2]+neu_state[3]*neu_state[3]
            pt = neu_state[4]*neu_state[4]+neu_state[5]*neu_state[5]
            
            p_active = pe+pm+pt
            
            # new xsections
            # Note : converting from natural units to [cm^2].            
            sig_NC = xs.nuDISxsection_NC_NusigmaInt(E_GeV,param)/(param.cm**2)
            sig_CC_e = xs.nuDISxsection_CC_NusigmaInt(E_GeV,0,param)/(param.cm**2)
            sig_CC_m = xs.nuDISxsection_CC_NusigmaInt(E_GeV,1,param)/(param.cm**2)
            sig_CC_t = xs.nuDISxsection_CC_NusigmaInt(E_GeV,2,param)/(param.cm**2)
            
            sig_total_active = sig_NC*p_active + sig_CC_e*pe + sig_CC_m*pm +sig_CC_t*pt
            
            p_survival = p_survival*np.exp(-pfac*sig_total_active)
            
            survival_prob_list.append([ineu,neu,E_GeV,r/(body.Radius*param.km),p_survival])
     
    if r >= track.xend or k>0.0:
        # if it is out or it had a CC-mu-e interaction
        if complete_decoherence :
            PMNS = MixMatrix(param)
            if not OPTNEUOSC_ERROR:
                result = ono.CompleteDecoherenceOscProb(ono.cvec(neu_state,param),PMNS,param)
            else :
                result = CompleteDecoherenceOscProb(cvec(neu_state,param),PMNS,param)    
        else : 
            if not OPTNEUOSC_ERROR:
                result = map(float,OscProb(ono.cvec(neu_state,param),param))
            else :
                result = map(float,OscProb(cvec(neu_state,param),param))
    else :
        # discart event
        result = []
        
    if not only_osc_prob:
        result.insert(0,E_GeV)
        result.insert(0,E_in)
        if neu == 1:
            result.insert(0,1)
        else:
            result.insert(0,0)    
        result.insert(0,ineu)
        result.insert(0,family)
        
        # adding also the survival probability.
        result.insert(5+param.numneu,p_survival)
    
    if file == None:
        if survival_prob:
            return survival_prob_list
        else :
            return result
    elif survival_prob:
        gt.hwritefile(file,survival_prob_list,param,header = False)
    else :
        gt.hwritefile(file,result,param,header = False)
        
    # cleaning up #
    if neu == 0 :
        param.neutype = "neutrino"
    elif neu == 1 :
        param.neutype = "antineutrino"
        
    del params        
    
#===============================================================================#
# INTERPOLATOR                                                                  #
#===============================================================================#
    
def InterOscProb(ineu,fneu,E,param,datapath = "../data/",Emin = 1.0,Emax = 1000000.0,return_interpolator = False,filename_append = ''):
    """ Returns P(nu_ineu -> nu_fneu) for osc. prob in the multi-GeV range from previously generated data by linear interpolation.
    
    Note : Uses data calculated by #L{DataNeuOscProb}. Updated to return the probability in any previouly
    calculated configuration.
    
    @type   ineu    :      integer
    @param  ineu    :      flavor of the initial neutrino (0 : e, 1 : mu, 2 : tau)
    @type   fneu    :      integer
    @param  fneu    :      flavor of the final neutrino (0 : e, 1 : mu, 2 : tau)
    @type   E       :      float
    @param  E       :      neutrino energy [eV]
    @type   param   :      PhysicsConstants
    @param  param   :      set of physical parameters to be used.     
    @type   datapath:      string
    @param  datapath:      path to the data files.
    @type   Emin    :      float
    @param  Emin    :      Minimum neutrino energy [GeV]
    @type   Emax    :      float
    @param  Emax    :      Maximum neutrino energy [GeV].
    @type   return_interpolator   :  Boolean
    @param  return_interpolator   :  If return interpolator = True : will return the interpolator instead of value.
    
    @rtype          :   float/interpolator
    @return         :   oscillation probability/interpolator
    """
    # constructing filename
    if filename_append == '':
        filename = "DataNeuOscProb_RK_"+param.neutype+"_Emin_"+str(Emin)+"_GeV_Emax_"+str(Emax)+"_GeV_ineu_"+str(ineu)+"_param_"+param.name+".dat"
    else : 
        filename = "DataNeuOscProb_RK_"+param.neutype+"_Emin_"+str(Emin)+"_GeV_Emax_"+str(Emax)+"_GeV_ineu_"+str(ineu)+"_param_"+param.name+"_"+str(filename_append)+".dat"
    
    # opening file and reading data
    file = open(datapath+filename,'r')
    h,dat = gt.hreadfilev3(file)
    dat = dat[0]
    # setuping up interpolator
    E_array = map(lambda x : x[0],dat)
    prob_array = map(lambda x : x[fneu+1],dat)
    inter = interpolate.interp1d(E_array,prob_array)
    # return interpolated value
    if return_interpolator:
        return inter
    else :
        return inter(E)
        
def InterSurProb(E,param,Emin = 1.0,Emax = 10000.0,datapath = "../data/NeutrinoAbsorptionProb/", return_interpolator = False):
    """ Returns the survival probability prob in the multi-GeV range from previously generated data by linear interpolation.

    @type   E       :       float
    @param  E       :       neutrino energy [eV]
    @type   param   :       physicsconstants
    @param  param   :       set of physical parameters to be used.
    @type   Emin    :       float
    @param  Emin    :       Minimum neutrino energy [GeV]
    @type   Emax    :       float
    @param  Emax    :       Maximum neutrino energy [GeV].
    @type   datapath:       string
    @param  datapath:       path to the data files.    
    @type   return_interpolator   :  Boolean
    @param  return_interpolator   :  If return interpolator = True : will return the interpolator instead of value.    
    """
    neutype = param.neutype
    name    = param.name
    
    if PC.pabs_n_inter == 0.0 or PC.pabs_a_inter == 0.0 :
        for neu in ["neutrino","antineutrino"]:
            filename = "DatNeuSunAbsorptionProbability_"+str(Emin)+"_GeV_"+str(Emax)+"_GeV_"+neu+".dat"
            
            file = open(datapath+filename,'r')            
            pabs = []
            gt.hreadfilev4(file,pabs,param)
            file.close()
            
            pabs = pabs[0]
            
            E_nu  = [p[0] for p in pabs]
            p_abs = [p[1] for p in pabs]
            
            inter_pabs = interpolate.interp1d(E_nu,p_abs)
            
            if neu == "neutrino":
                PC.pabs_n_inter = inter_pabs
            elif neu == "antineutrino":
                PC.pabs_a_inter = inter_pabs
                
    param.neutype = neutype
    param.name = name
    if param.neutype == "neutrino":
        inter_pabs = PC.pabs_n_inter
    elif param.neutype == "antineutrino":
        inter_pabs = PC.pabs_a_inter
    else:
        print "Wrong neutype."
        quit()
    
    if return_interpolator:
        return inter_pabs
    else :
        return inter_pabs(E)            

#===============================================================================
# DATA GENERATORS
#===============================================================================

def GenDataOscProbFlavor(ineu,Emin,Emax,param,body = bd.Vacuum(),track = None,savepath = "../data/test/",binnum = 100, integrator = "GSL",complete_decoherence = False,filename_append = ''):
    """ Calculates and saves the probabilities for a given configuration file for a single neutrino flavor
    specified by ineu.

    @type   Emin    :   float
    @param  Emin    :   minimum neutrino energy [eV]
    @type   Emax    :   float
    @param  Emax    :   maximum neutrino energy [eV]
    @type   body    :   body
    @param  body    :   body with the asociated density profile.
    @type   param   :   physicsconstants
    @param  param   :   set of physical parameters to be used.
    @type   savepath:   string
    @param  savepath:   location where the data will be saved    
    @type   binnum  :   integer
    @param  binnum  :   number of bins
    @type   integrator  :   string
    @param  integrator  :   "GSL": GSL integrator, "Py" : python integrator
    @type   complete_decoherence    :   boolean
    @param  complete_decoherence    :   toogle complete decoherence ON/OFF
    @type   filename_append    :   string
    @param  filename_append    :   append string to file name.
    """
    param.Refresh()
    fM2 = flavorM2(param)
    PMNS = MixMatrix(param)
    
    if track == None :
        Ri = 0.001*body.Radius*param.km
        Rf = body.Radius*param.km
        track = body.track(Ri,Ri,Rf)
    
    ERK = gt.LogSpaceEnergies(Emin,Emax,binnum = binnum)

    if integrator == "GSL":
        oscprob = map(lambda E : CalNeuOscGSL(ineu,E,track,body,fM2,param,file = None,family = 0,NCCC_int = False,complete_decoherence = complete_decoherence),ERK)
    elif integrator == "Py":
        oscprob = map(lambda E : CalNeuOscPy(ineu,E,track,body,fM2,param,complete_decoherence = complete_decoherence),ERK)
    else :
        print "NC:NEU:NEUOSC:ERROR: GenDataOscProbFlavor : wrong integrator."
        quit()
    
    pp  = []
    delta = max(len(oscprob[0])-param.numneu,0)
    
    for i,p in enumerate(oscprob):
        tmp = []
        tmp.append(ERK[i])
        for j in range(param.numneu):
            tmp.append(p[delta+j])
        pp.append(tmp)
        
    if filename_append == '':
        filename = "DataNeuOscProb_RK_"+param.neutype+"_Emin_"+str(Emin/param.GeV)+"_GeV_Emax_"+str(Emax/param.GeV)+"_GeV_ineu_"+str(ineu)+"_param_"+param.name+".dat"
    else : 
        filename = "DataNeuOscProb_RK_"+param.neutype+"_Emin_"+str(Emin/param.GeV)+"_GeV_Emax_"+str(Emax/param.GeV)+"_GeV_ineu_"+str(ineu)+"_param_"+param.name+"_"+str(filename_append)+".dat"
        
    file = open(savepath+filename,'w')
    gt.hwritefile(file,pp,param)
    file.close()
    
    
def GenDataOscProb(Emin, Emax, param, binnum = 100, zip = False, body = bd.Vacuum(), track = None, savepath = "../data/", integrator = "GSL",complete_decoherence = False,filename_append = ''):
    """ Calculates and saves the probabilities for a given configuration file and all initial neutrino
    flavors. For the same calculation for a given neutrino flavor look for : GenDataOscProbFlavor.

    @type   Emin    :   float
    @param  Emin    :   minimum energy [GeV]
    @type   Emax    :   float
    @param  Emax    :   maximum energy [GeV]
    @type   param   :   physicsconstants
    @param  param   :   set of physical parameters to be used.
    @type   binnum  :   integer
    @param  binnum  :   number of bins
    @type   zip     :   boolean
    @param  zip     :   if true will zip all files together
    @type   body    :   body
    @param  body    :   body where the neutrino osc. will take place.
    @type   track   :   body.track
    @param  track   :   trayectory of the neutrino in the body.
    @type   savepath:   string
    @param  savepath:   location where the data will be saved
    @type   integrator  :   string
    @param  integrator  :   "GSL": GSL integrator, "Py" : python integrator
    @type   complete_decoherence    :   boolean
    @param  complete_decoherence    :   toogle complete decoherence ON/OFF    
    @type   filename_append    :   string
    @param  filename_append    :   append string to file name.    
    """
    Emin = Emin*param.GeV
    Emax = Emax*param.GeV
    
    # generate neutrino/antineutrino probabilities
    
    print "======= BEGIN "+body.name+" OSC PROB CALC USING RK ======="
    
    param.neutype = "neutrino"
    param.Refresh()
    for ineu in [0,1,2]:
        print "BEGIN : CALCULATION : PROB : CONF = "+param.name+" : INEU = "+str(ineu)+" : NEUTYPE = "+param.neutype
        GenDataOscProbFlavor(ineu,Emin,Emax,param,
                             body = body,track = track, savepath = savepath,binnum = binnum,integrator = integrator,
                             complete_decoherence = complete_decoherence, filename_append = filename_append)
        print "END : CALCULATION : PROB : CONF = "+param.name+" : INEU = "+str(ineu)
        
    param.neutype = "antineutrino"
    param.Refresh()
    for ineu in [0,1,2]:
        print "BEGIN : CALCULATION : PROB : CONF = "+param.name+" : INEU = "+str(ineu)+" : NEUTYPE = "+param.neutype
        GenDataOscProbFlavor(ineu,Emin,Emax,param,
                             body = body, track = track, savepath = savepath,binnum = binnum,integrator = integrator,
                             complete_decoherence = complete_decoherence, filename_append = filename_append)
        print "END : CALCULATION : PROB : CONF = "+param.name+" : INEU = "+str(ineu)
        
    print "Conf. "+ param.name +" probabilities generated."
    
    print "======= END SUN OSC PROB CALC USING RK ======="
    
    if zip :
        print "======= BEGIN ZIP CREATION ======="
        import zipfile as zp
        ## will put files in zip file ##
        zipfilename = "output.zip"
        zipfile = zp.ZipFile(savepath+zipfilename, mode = 'a')
        
        filenames = []
        for ff in os.listdir(savepath):
            if fnmatch.fnmatch(ff,"*"+param.name+"*.dat"):
                filenames.append(ff)
                
        for filename in filenames :
            zipfile.write(savepath+filename,filename)
            
        zipfile.close
        
        print "======= END ZIP CREATION ======="
        
def GenDataSurProb(Emin,Emax,body = bd.Sun(),param = PC.PhysicsConstants()):
    """ Calculates and saves the neutrino survival probability.
    
    UPDATED TO INCLUDE CC+NC INTERACTIONS IN ABSORPTION.

    @type   Emin    :   float
    @param  Emin    :   neutrino energy [eV]
    @type   Emax    :   float
    @param  Emax    :   neutrino energy [eV]    
    @type   body    :   body
    @param  body    :   body where the neutrino will propagate.
    @type   param   :   physicsconstants
    @param  param   :   set of physical parameters to be used.
    """
    import body as bd
    Sun   = bd.Sun()
    R_sun = Sun.Radius
    
    if param.neutype == "neutrino":
        neu = 0
    elif param.neutype == "antineutrino":
        neu = 1
    else:
        print "NC:NEU:NEUOSC:ERROR: GenDataSurProb : wrong integrator."
        quit()
        
    E_range = gt.LogSpaceEnergies(Emin,Emax, binnum = 200)
    p_abs   = []    
    
    for E in E_range:
        integrant = integrate.quad(lambda r :
                        #(xs.signuNCC(E/param.GeV,neu))*
                        (xs.signuNCC(E/param.GeV,neu)+xs.signuNNC(E/param.GeV,neu))*                        
                        (Sun.rdensity(r/R_sun)*(param.gr/param.cm**3)/(param.atomic_mass_unit*param.gr))*
                        (param.km),0.0,R_sun)[0]
        
        p_abs.append([E,
                      np.exp(-integrant)
                     ])
    
    datapath = "../data/NeutrinoAbsorptionProb/"
    filename = "DatNeuSunAbsorptionProbability_"+str(Emin/param.GeV)+"_GeV_"+str(Emax/param.GeV)+"_GeV_"+str(param.neutype)+".dat"
    file = open(datapath+filename,'w')
    
    gt.hwritefile(file,p_abs,param)
    
    print "Files created."
    
def GenEvents(ineu,Emin,Emax,num_events,body = bd.Sun(),param= PC.PhysicsConstants(),complete_decoherence = False):
    """ Generates events using #L{CalNeuOscGSL}.

    @type   ineu    :   int
    @param  ineu    :   initial neutrino state
    @type   Emin    :   float
    @param  Emix    :   neutrino energy [GeV]
    @type   Eman    :   float
    @param  Emax    :   neutrino energy [GeV]    
    @type   num_events   :   int
    @param  num_events   :   number of events generated
    @type   param   :   physicsconstants
    @param  param   :   set of physical parameters to be used.
    """
    datapath = param.savedatapath
        
    PMNS    = MixMatrix(param)
    fM2     = flavorM2(param)

    # file handling
    if param.savefilename == None:
        now = datetime.datetime.now()
     
        n = 1
        for ff in os.listdir(datapath):
            if fnmatch.fnmatch(ff,'evt_MC_ineu_'+str(ineu)+"_"+param.name+"_"+str(now.day)+"_"+str(now.month)+"_"+str(now.year)+"*.dat"):
		n = n + 1

            if n == 1 :
                filename    = "evt_MC_ineu_"+str(ineu)+"_"+param.name+"_"+str(now.day)+"_"+str(now.month)+"_"+str(now.year)+".dat"
            else :
                filename    = "evt_MC_ineu_"+str(ineu)+"_"+param.name+"_"+str(now.day)+"_"+str(now.month)+"_"+str(now.year)+"."+str(n)+".dat"
 
    else: 
        filename = param.savefilename
    
    file        = open(datapath+filename,'ar')
    
    gt.hwritefile(file,[],param,header = True)

    # build log binned energy list
    E_list = [Emin]
    if Emin == Emax :
        E_list = [Emin]
    else :
        E_list = gt.MidPoint(gt.LogSpaceEnergies(Emin,Emax))
    
    # starting up family number
    family      = 1
    
    # generating events
    for EE in E_list:
        if EE>=5.0 and EE<=1000.0:
            for i in range(num_events):
                Ri      = 0.01*Sun.Radius*param.km
                Rf      = 1.0*Sun.Radius*param.km
                track   = Sun.track(Ri,Ri,Rf)
                CalNeuOscGSL(ineu,EE*param.GeV,track,body,fM2,PMNS,param,file,family,True,complete_decoherence = complete_decoherence)
                family = family + 1      
                
    print "NEUOSC:DONE: Finish generating events."  
        
#===============================================================================
# DENSITY MATRIX FORMALISM - ON DEVELOPMENT - NOT FULLY TESTED
#===============================================================================

def flavorS(L,track,param,flavorM2,body,E):
    """ Calculates evolution operator in constant density and molarity
    using Cayley-Hamilton formalism - omitting a global phase factor.

    @type   L       :   float
    @param  L       :   baseline [eV^-1]
    @type   track   :   physicsconstants
    @param  track   :   set of physical parameters to be used.
    @type   param   :   string
    @param  param   :   location where the data will be saved
    @type   flavorM2:   boolean
    @param  flavorM2:   if true will zip all files together    
    @type   body    :   body
    @param  body    :   body where the neutrino propagates.
    @type   E       :   float
    @param  E       :   neutrino energy [eV]
    """
    # L : propagating distance in eV^-1
    S = np.zeros([param.numneu,param.numneu],complex)
    #print S
    ii = complex(0.0,-1.0)
    fH = flavorH(track,param,flavorM2,body,E)/ii
    Id = np.identity(param.numneu,complex) 
    T = fH - np.trace(fH)*Id/3.0
    print np.trace(T)
    c1 = - np.trace(np.dot(T,T))/2.0
    eigen,V = gt.eigenvectors(T)
    reigen = map(lambda x : x.real, eigen)
    ieigen = map(lambda x : x.imag, eigen)
    
    for i in np.arange(0,param.numneu,1):
    	print (np.exp(ii*eigen[i]*L)/(3*eigen[i]**2+c1))
        S = S + (np.exp(ii*eigen[i]*L)/(3*eigen[i]**2+c1))*((eigen[i]**2+c1)*Id+eigen[i]*Id+np.dot(T,T))
    return S

def NuEvolution(ineu,fneu,E,track,body,flavorM2,param):
    """ Evolves a state a number of steps.
    
    Note : Incomplete.

    @type   ineu    :   integer 
    @param  ineu    :   initial neutrino flavor
    @type   fneu    :   integer 
    @param  fneu    :   final neutrino flavor
    @type   E       :   float
    @param  E       :   neutrino energy [eV]    
    @type   track   :   physicsconstants
    @param  track   :   set of physical parameters to be used.
    @type   body    :   body
    @param  body    :   body where the neutrino propagates.    
    @type   flavorM2:   boolean
    @param  flavorM2:   if true will zip all files together    
    @type   param   :   string
    @param  param   :   location where the data will be saved
    """
    Ri = track.Ri
    Rf = track.Rf
    nstep = 10.0
    step = (Rf - Ri)/nstep
    S = np.identity(param.numneu, complex)
    for i in np.arange(0,nstep,1):
        track.x = track.x + step
        S = np.dot(flavorS(step,track,param,flavorM2,body,E),S)
    return S        

#===============================================================================
# OTHER TOOLS
#===============================================================================

def AdiabaticProbability(ineu,fneu,E,Ri,body,PMNS,fM2,param):
    Rf = body.Radius*param.km
    track = body.track(Ri,Ri,Rf)    
    ##
    ii = complex(0.0,-1.0)
    fH = flavorH(track,param,fM2,body,E)/ii
    ##
    D, V = gt.eigenvectors(fH)
    vfinal = PMNS.U[fneu,:]
    prob = 0.0
    for i in np.arange(0,param.numneu,1):
            prob = prob + gt.norm(V[i][fneu])**2*gt.norm(vfinal[i])**2
    return prob

def Eeigenvals(E,body,track,fM2,param):
    """ Calculates the Eigenvalues of the hamiltonian including the mass term.
    This requieres a body and a position in the body given by a body.track object.
    
    @type   E       :   float
    @param  E       :   neutrino energy [eV]
    @type   body    :   body
    @param  body    :   body
    @type   track   :   body.track object
    @param  track   :   location inside the body
    @type   fM2     :   matrix
    @param  fM2     :   mass matrix in the flavor basis
    @type   param   :   physicsconstants
    @param  param   :   set of physics parameters to be used.
    """
    param.Refresh()
    ##
    ii = complex(0.0,-1.0)
    fH = flavorH(track,param,fM2,body,E)/ii
    ##
    D = np.linalg.eigvals(fH)
    DR = map(lambda d : d.real,D)
    return sorted(DR)

def NeuComposition(i,E,x,body,fM2,param,flavor):
    """ Returns composition of neutrino i = 1,2,3,...,numneu (e,mu,tau,s1,s2...,s_numeu-3)
    
    Note : i starts from 0.

    @type   i       :   float
    @param  i       :   minimum energy [GeV]
    @type   E       :   float
    @param  E       :   maximum energy [GeV]
    @type   x       :   physicsconstants
    @param  x       :   set of physical parameters to be used.
    @type   body    :   integer
    @param  body    :   number of bins
    @type   fM2     :   boolean
    @param  fM2     :   if true will zip all files together
    @type   param   :   string
    @param  param   :   location where the data will be saved
    @type   flavor  :   boolean
    @param  flavor  :   if flavor = True : returns composition as a function of the mass basis.
                        if flavor = False: returns composition as a function of the flavor basis.
    """
    # starting up mixing matrix
    PMNS = MixMatrix(param)
    U = PMNS.U
    UCT = PMNS.UCT
    # defining a body track
    Ri = x*body.Radius*param.km
    Rf = body.Radius*param.km
    track = body.track(Ri,Rf)    
    #
    ii = complex(0.0,-1.0)
    fH = flavorH(track,param,fM2,body,E)/ii
    ##
    D, V = gt.eigenvectors(fH)
    if flavor :
        vmf = []
        for j in np.arange(0,param.numneu,1):
            vmf.append(map(lambda x : gt.norm(x)**2, V[:][j]))
        vv = map(lambda x : x[i] ,vmf)
        return vv
    else : 
        vv = map(lambda x : gt.norm(x)**2, V[:][i])
        return vv      
        
if __name__ == '__main__':
    pass
