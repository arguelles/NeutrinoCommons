""" Small library with the optimiced python osc. calculation functions.
This library is implemented in Cython. And is meant to be fully compatible
with #{neuosc}.
"""    
from libc.math cimport sin,cos,sqrt,log
from numpy import zeros,dot,array
cimport numpy as np

from scipy.linalg.matfuncs import expm2

# importing some UNoptimice libraries

import neutrinocommon.physconst.physicsconstants as PC
import neutrinocommon.neu.xsections as xs
import neutrinocommon.neu.neuosc as no

# importing optimzed libraries
try :
    import neutrinocommon.shareobjects.optgeneraltools as gt
except : 
    import neutrinocommon.tools.generaltools as gt
try : 
    import neutrinocommon.shareobjects.taudecay as td
except :
    pass
    
# python modules
import numpy as np
import pygsl.odeiv as odeiv
import scipy as sp

# FUNCTIONS #

cdef extern from "complex.h":
    double norm(complex)

cdef double norm(complex x):
    return sqrt(x.real*x.real + x.imag*x.imag)
    
cdef double chi2(double th, double exp):
    cdef double chi
    chi = (exp - th)*(exp - th)/exp
    return chi
    
cpdef cvec(vv,param):
    cdef int i
    neu = zeros([param.numneu,1],complex)
    for i in range(0,2*param.numneu,2):
        neu[i/2,0] = complex(vv[i],vv[i+1])
    return neu

cpdef rvec(vv,param):
    cdef int i
    neu = np.zeros([2*param.numneu],float)
    for i in np.arange(0,2*param.numneu,2):
        neu[i] = vv[i/2,0].real
        neu[i+1] = vv[i/2,0].imag
    return neu
    
cpdef OscProb(cvec,param):
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

cpdef flavorAcc(param,E,body,track):
    # Charge-current potential
    # E : neutrino energy
    # body : body where the neutrino propagetes
    # track : trayectory
    Acc = zeros([param.numneu,param.numneu],complex)
    ye = body.ye(track)
    CC = param.sqr2*param.GF*param.Na*param.cm**-3*body.density(track)*ye
    #NC = CC*body.RNC(track)
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
        
cpdef flavorH(track,param,flavorM2,body,E):
    # Flavor Hamiltonian WARNING FACTOR II in front
    # rhs of the Mikheev-Smirnov-Wolfenstein (MSW) eq.
    ii = complex(0.0,-1.0)
    return (ii)*(flavorM2/(2.0*E)+flavorAcc(param,E,body,track))

cpdef RHS(neuvec,double x,track,param,flavorM2,body,double E):
    # Right hand side in Schrodinger's equations
    track.x = x
    ## this should be done in complex plane -- TMP
    neuvec = cvec(neuvec,param)
    rhs = dot(flavorH(track,param,flavorM2,body,E),neuvec)
    tmprhs = (array([[val[0].real,val[0].imag] for val in rhs])).flatten()        
    return tmprhs

cpdef RHS_GSL(double x,neuvec,params):
    track,param,flavorM2,body,E = params.track, params.param,params.flavorM2,params.body,params.E
    # Right hand side in Schrodinger's equations
    track.x = x
    ## this should be done in complex plane -- TMP
    neuvec = cvec(neuvec,param)
    rhs = dot(flavorH(track,param,flavorM2,body,E),neuvec)
    tmprhs = (array([[val[0].real,val[0].imag] for val in rhs])).flatten()        
    return tmprhs

cpdef R(int i, int j,int cp,param):
    cdef int k
    cdef double sd,cd
    cdef complex faseCP
        
    # rotation matrix generator
    R = zeros([param.numneu,param.numneu],complex)
    # diagonal terms
    for k in range(param.numneu):
        if(k != i-1 and k != j-1):
            R[k,k] = 1.0
        else :
            R[k,k] = param.c[i,j]
    # non-diagonal terms
    if(cp != 0):
        sd = sin(param.dcp[cp])
        cd = cos(param.dcp[cp])
        faseCP = complex(cd,sd)
    else:
        faseCP = complex(1.0,0.0)
        
    R[i-1,j-1] = param.s[i,j]*faseCP.conjugate()
    R[j-1,i-1] = -param.s[i,j]*faseCP
    return R
    
cdef Ropt(int i,int j,int nn,int cp,double s[4][4],double c[4][4],double dcp[4]):
    cdef int k
    cdef double sd,cd
    cdef complex faseCP
        
    # rotation matrix generator
    R = zeros([nn,nn],complex)
    # diagonal terms
    for k in range(nn):
        if(k != i-1 and k != j-1):
            R[k,k] = 1.0
        else :
            R[k,k] = c[i][j]
    # non-diagonal terms
    if(cp != 0):
        sd = sin(dcp[cp])
        cd = cos(dcp[cp])
        faseCP = complex(cd,sd)
    else:
        faseCP = complex(1.0,0.0)
        
    R[i-1,j-1]  == s[i][j]*faseCP.conjugate()
    R[j-1,i-1] = -s[i][j]*faseCP
    return R    

cpdef CompleteDecoherenceOscProb(fneu_flavor,PMNS,param):
    cdef int flv
    fneu_mass = dot(PMNS.UCT,fneu_flavor)
    prob_flv = []
    for flv in range(param.numneu):
        vfinal = PMNS.U[flv,:]
        prob = 0.0
        for i in range(param.numneu):
            prob = prob + gt.norm(fneu_mass[i])**2*gt.norm(vfinal[i])**2
        prob_flv.append(prob[0])
    return prob_flv


#===============================================================================
# FAST DATA COUNTING
#===============================================================================
    
cpdef DMEvtProcessv2(E_nu_bin,E_anu_bin,E_nu_list,E_bin_width,events,DM_pdf_table,double E_nu_list_0,double log_E_bin_ratio,double DMm_GeV, double E_in_GeV,int ini,fam_num,ineu):
    cdef int i,ii,jj,iineu,ineutype,family
    cdef double E_nu_in,E_nu_out
    
    cdef double eps
    
    eps = 1.0e-4
    
    for i,e in enumerate(events):
        
        if len(e) > 5:
            family = e[0]
            iineu  = e[1]
            ineutype = e[2]
            E_nu_in  = e[3]
            E_nu_out = e[4]
            
            fam_factor = fam_num[family-1]
            
            ini_evt_energy = E_in_GeV
            
            ii = int(log(E_nu_out/E_nu_list_0)/log_E_bin_ratio+eps)
            
            #if E_nu_out != E_in_GeV :
            #    print e
            #    print E_nu_out, e[6], ii
            #    raw_input("Wait ...")
            #ini = int(log(ini_evt_energy/E_nu_list_0)/log_E_bin_ratio+eps)
            
            iineu = ineu
            
            if ineutype == 0:
                DM_pdf_array = DM_pdf_table[iineu*2]
                if E_nu_in + E_bin_width[ini] > DMm_GeV:
                    E_nu_bin[ii]  = E_nu_bin[ii]  + e[6]*(DM_pdf_array[ini])*(DMm_GeV-E_nu_list[ini])/fam_factor
                else:
                    E_nu_bin[ii]  = E_nu_bin[ii]  + e[6]*(DM_pdf_array[ini])*E_bin_width[ini]/fam_factor
            elif ineutype == 1 :
                DM_pdf_array = DM_pdf_table[iineu*2+1]
                if E_nu_in + E_bin_width[ini] > DMm_GeV :
                    E_anu_bin[ii] = E_anu_bin[ii] + e[6]*(DM_pdf_array[ini])*(DMm_GeV-E_nu_list[ini])/fam_factor
                else :
                    E_anu_bin[ii] = E_anu_bin[ii] + e[6]*(DM_pdf_array[ini])*E_bin_width[ini]/fam_factor
            
    return [E_nu_bin,E_anu_bin]    
    
cpdef DMEvtProcess_Fluxv2(E_nu_bin,E_anu_bin,E_nu_list,E_bin_width,events,DM_pdf_table,double E_nu_list_0,double log_E_bin_ratio,double DMm_GeV,E_in_GeV,fam_num,ineu):    
    cdef int i,ii,jj,ini,iineu,ineutype,family
    cdef double E_nu_in,E_nu_out
    
    cdef double eps
    
    eps = 1.0e-4
    
    for i,e in enumerate(events):
        
        if len(e) > 5:
            family = e[0]
            iineu = e[1]
            ineutype = e[2]
            E_nu_in  = e[3]
            E_nu_out = e[4]
            
            fam_factor = fam_num[family-1]
            ini_evt_energy = E_in_GeV
            
            iineu = ineu
            
            ii = int(log(E_nu_out/E_nu_list_0)/log_E_bin_ratio+eps)
            ini = int(log(ini_evt_energy/E_nu_list_0)/log_E_bin_ratio+eps)
            
            if ini != ii :
                ii = ini
                if ineutype == 0:
                    DM_pdf_array = DM_pdf_table[iineu*2]
                    E_nu_bin[ii]  = E_nu_bin[ii]  + e[6]*(DM_pdf_array[ini])/fam_factor
                elif ineutype == 1 :
                    DM_pdf_array = DM_pdf_table[iineu*2+1]
                    E_anu_bin[ii] = E_anu_bin[ii] + e[6]*(DM_pdf_array[ini])/fam_factor
            
    return [E_nu_bin,E_anu_bin]        

#===============================================================================
# Neu Osc Rediagonalizing Calculation - Cython
#===============================================================================

cpdef RHS_INT_GSL(double x,neuvec,params):
    cdef complex ii
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

    interaction_H = dot(S,dot(H1,S.conjugate().T))
    
    rhs_int = (-ii)*interaction_H
    
    rhs = dot(rhs_int,neuvec)    
    tmprhs = (array([[val[0].real,val[0].imag] for val in rhs])).flatten()
    
    #storing stuff
    params.cH = H_current
    params.cS = S
    return tmprhs

cdef class ParamsContainer:
    cdef object track,param,flavorM2,body,E,H0,S,rpos,cH,cS

    property track:
        def __get__(self):
            return self.track
        def __set__(self,value):
            self.track = value
        def __del__(self):
            track = None
            
    property param:
        def __get__(self):
            return self.param
        def __set__(self,value):
            self.param = value
        def __del__(self):
            param = None            

    property flavorM2:
        def __get__(self):
            return self.flavorM2
        def __set__(self,value):
            self.flavorM2 = value
        def __del__(self):
            flavorM2 = None            
            
    property body:
        def __get__(self):
            return self.body
        def __set__(self,value):
            self.body = value
        def __del__(self):
            body = None            
            
    property E:
        def __get__(self):
            return self.E
        def __set__(self,value):
            self.E = value
        def __del__(self):
            E = None            
            
    property H0:
        def __get__(self):
            return self.H0
        def __set__(self,value):
            self.H0 = value
        def __del__(self):
            H0 = None            

    property S:
        def __get__(self):
            return self.S
        def __set__(self,value):
            self.S = value
        def __del__(self):
            S = None                        
            
    property rpos:
        def __get__(self):
            return self.rpos
        def __set__(self,value):
            self.rpos = value
        def __del__(self):
            rpos = None                        
            
    property cH:
        def __get__(self):
            return self.cH
        def __set__(self,value):
            self.cH = value
        def __del__(self):
            cH = None                        

    property cS:
        def __get__(self):
            return self.cS
        def __set__(self,value):
            self.cS = value
        def __del__(self):
            cS = None                        

cpdef CalNeuOscGSL(ineu,E,track,body,flavorM2,param,file = None,family = 0,NCCC_int = False,complete_decoherence = False,survival_prob = False,E_threshold = 1.0,optimization = True,only_osc_prob = True,abs_error = 1.0e-7, rel_error = 1.0e-6):
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
    cdef complex ii
    cdef int dim,m,n,k
    cdef double h,r,E_in,p_survival,r_old,r_new,E_MC,p_MC
    
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
    params.H0 = ii*flavorH(track,param,flavorM2,body,E)
    params.S  = np.identity(param.numneu, complex)
    
    # setup GSL solver
    dim     = int(2*param.numneu)
    stepper = odeiv.step_rkf45
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
            
            dsde_NC     = xs.dsde_NC_opt(E_GeV,E_MC,param)
            dsde_CC_e   = xs.dsde_CC_opt(E_GeV,E_MC,param,flv = 0)
            # Assuming that the muon diff. cross section is similar to the electron
            dsde_CC_m   = dsde_CC_e
            # else we should use this line :
            #dsde_CC_mu = xs.dsde_CC_opt(E_GeV,E_MC,param,flv = 1)
            dsde_CC_t   = xs.dsde_CC_opt(E_GeV,E_MC,param,flv = 2)
            
            # new xsections
            sig_NC = xs.nuDISxsection_NC_NusigmaInt(E_GeV,param)
            sig_CC_e = xs.nuDISxsection_CC_NusigmaInt(E_GeV,0,param)
            sig_CC_m = xs.nuDISxsection_CC_NusigmaInt(E_GeV,1,param)
            sig_CC_t = xs.nuDISxsection_CC_NusigmaInt(E_GeV,2,param)
            
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
        else:
            E_GeV = params.E/param.GeV
            # calculate probability survival
            pfac = np.abs((body.density(track)/param.atomic_mass_unit)*(r_adv/param.cm))
            
            pe = neu_state[0]*neu_state[0]+neu_state[1]*neu_state[1]
            pm = neu_state[2]*neu_state[2]+neu_state[3]*neu_state[3]
            pt = neu_state[4]*neu_state[4]+neu_state[5]*neu_state[5]
            
            p_active = pe+pm+pt
            
            # new xsections
            sig_NC = xs.nuDISxsection_NC_NusigmaInt(E_GeV,param)
            sig_CC_e = xs.nuDISxsection_CC_NusigmaInt(E_GeV,0,param)
            sig_CC_m = xs.nuDISxsection_CC_NusigmaInt(E_GeV,1,param)
            sig_CC_t = xs.nuDISxsection_CC_NusigmaInt(E_GeV,2,param)
            
            sig_total_active = sig_NC*p_active + sig_CC_e*pe + sig_CC_m*pm +sig_CC_t*pt
            
            p_survival = p_survival*np.exp(-pfac*sig_total_active)
            
            if survival_prob :
                survival_prob_list.append([ineu,neu,E_GeV,r/(body.Radius*param.km),p_survival])
     
    if r >= track.xend or k>0.0:
        # if it is out or it had a CC-mu-e interaction
        if complete_decoherence :
            PMNS = no.MixMatrix(param)
            result = CompleteDecoherenceOscProb(cvec(neu_state,param),PMNS,param)    
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
        
    params = ParamsContainer()
