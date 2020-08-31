"""
Author  : C.A. Arguelles
Date    : 10/MAY/2011

This script contains information of the neutrino
cross sections, including differential cross sections.
It also implements the NUSIGMA interface.

Log:
- Modified on 9/FEB/2012 : adding routines by M. Bustamente
to include neutrion-nucleus interactions.

- Modified on 15/ABR/2012 by C.Arguelles:
    + All functions will now return the cross section, or differential
    cross section, in natural units.

"""

import scipy.interpolate as interpolate
import numpy as np
import os
import subprocess
# my module
import neutrinocommon
import neutrinocommon.physconst.physicsconstants as PC
import neutrinocommon.tools.generaltools as gt
try : 
    import neutrinocommon.shareobjects.nudsde as oxs
    NUDSDE_ERROR = False
except ImportError:
    print "NC:NEU:XSECTIONS:ERROR: Loading NUSIGMA interface : nudsde."
    NUDSDE_ERROR = True

# global variables
filepath = neutrinocommon.neu.__path__

global_datapath = filepath[0]+"/data/xs/"
global_binpath  = filepath[0]+"/bin/xs/"

pc = PC.PhysicsConstants()

#### THIS IS TEMPORARY ####
# Number density of deuterium atoms
nD = 2.2149e22 # [cm^-3]
# Number density of Fe atoms
nFe = 8.4903e22 # [cm^-3]
# Number density of Pb atoms
nPb = 3.2959e22 # [cm^-3]
# Number density of C atoms
nC = 6.0271599e23 # [cm^-3]

#===============================================================================#
# Differential Cross Sections                                                   #
#===============================================================================# 

#####################################
##    USING nusigma-FORTRAN bin    ##
#####################################

def dsde_CC_NusigmaBin(E_1,E_2,param):
    """ Uses the FORTRAN bin to calculate the neutrino/antineutrino
    differential charge current cross section. The FORTRAN executable
    was generated using the NUSIGMA FORTRAN library.
        
    Implements:
    
    ds/dE in nu + N -> l + X,
    
    assuming isoscalar target. 
    
    Note: This function calls DIRECTLY the FORTRAN bin. Function
    replaced by the faster NUSIGMA-interface.
    
    @type   E_1 	:	float
    @param  E_1	    :	neutrino energy [GeV]
    @type   E_2 	:	float
    @param  E_2	    :	lepton energy [GeV]
    @type   param   :	physics parameter set
    @param  param   :	physics parameter set

    @rtype	        :	float
    @return	        :	neutrino differential CC cross section. New : [eV^-3] Old : [cm^2 GeV^-1]
    """
    if param.neutype == "neutrino":
        return float(subprocess.Popen(global_binpath+'nudsde_sp.exe '+str(E_1)+ ' ' +str(E_2)+' ' +str(1)+' N '+' CC',shell=True,stdout=subprocess.PIPE).communicate()[0])*(param.cm**2/param.GeV)
    elif param.neutype == "antineutrino":
        return float(subprocess.Popen(global_binpath+'nudsde_sp.exe '+str(E_1)+ ' ' +str(E_2)+' ' +str(2)+' N '+' CC',shell=True,stdout=subprocess.PIPE).communicate()[0])*(param.cm**2/param.GeV)
    else:
        print "NC:NEU:XSECTIONS:ERROR: dsde_CC_NusigmaBin : Wrong neutrino type."
        quit()

def dsde_NC_NusigmaBin(E_1,E_2,param):
    """ Uses the FORTRAN bin to calculate the neutrino/antineutrino
    differential neutral current cross section. The FORTRAN executable
    was generated using the NUSIGMA FORTRAN library.
        
    Implements:
    
    ds/dE in nu + N -> nu + N,
    
    assuming isoscalar target. 
    
    Note: This function calls DIRECTLY the FORTRAN bin. Function
    replaced by the faster NUSIGMA-interface.
    
    @type   E_1 	:	float
    @param  E_1	    :	neutrino energy [GeV]
    @type   E_2 	:	float
    @param  E_2	    :	lepton energy [GeV]
    @type   param   :	physics parameter set
    @param  param   :	physics parameter set

    @rtype	        :	float
    @return	        :	neutrino differential NC cross section. New : [eV^-3] Old : [cm^2 GeV^-1]
    """
    if param.neutype == "neutrino":
        return float(subprocess.Popen(global_binpath+'nudsde_sp.exe '+str(E_1)+ ' ' +str(E_2)+' ' +str(1)+' N '+' NC',shell=True,stdout=subprocess.PIPE).communicate()[0])*(param.cm**2/param.GeV)
    elif param.neutype == "antineutrino":
        return float(subprocess.Popen(global_binpath+'nudsde_sp.exe '+str(E_1)+ ' ' +str(E_2)+' ' +str(2)+' N '+' NC',shell=True,stdout=subprocess.PIPE).communicate()[0])*(param.cm**2/param.GeV)
    else:
        print "NC:NEU:XSECTIONS:ERROR: dsde_NC_NusigmaBin : Wrong neutrino type."
        quit()        
    
    
def dsde_CC_NusigmaBin_ListOpt(E_1,E_2,param,flv = 0):
    """ Uses the FORTRAN bin to calculate the neutrino/antineutrino
    differential charge current cross section. The FORTRAN executable
    was generated using the NUSIGMA FORTRAN library. Optimized for
    several calls of the same E_1.
    
    ds/dE in nu + N -> l + X
    
    Assuming isoscalar target. The neutrino flavor is specified by flv.
    
    @type   E_1 	:	float
    @param  E_1 	:	neutrino energy [GeV]
    @type   E_2 	:	float
    @param  E_2 	:	lepton energy [GeV]
    @type   param   :	physics parameter set
    @param  param   :	physics parameter set
    @type   flv     :	integer
    @param  flv     :	0 : e, 1 : mu, 2 : tau

    @rtype	        :	float
    @return	        :	neutrino differential CC DIS cross section. New : [eV^-3] Old : [cm^2 GeV^-1]
    """
    if PC.act_dsde_CCe_inter == 0 or PC.act_dsde_CCm_inter == 0 or PC.act_dsde_CCt_inter == 0 or PC.E_CC_act != E_1 :
        for flavor in [0,1,2]:
            if param.neutype == "neutrino":
                datarray = subprocess.Popen(binpath+'nudsde_sp_list.exe '+str(E_1)+ ' ' +str(2*flavor+1)+' N '+' CC',shell=True,stdout=subprocess.PIPE).communicate()[0]*(param.cm**2/param.GeV)
            elif param.neutype == "antineutrino":
                datarray = subprocess.Popen(binpath+'nudsde_sp_list.exe '+str(E_1)+ ' ' +str(2*flavor+2)+' N '+' CC',shell=True,stdout=subprocess.PIPE).communicate()[0]*(param.cm**2/param.GeV)
            else:
                print "Wrong neutrino type"
                quit()
            datarray = datarray.split()
            E_lep = [float(datarray[i]) for i in range(0,len(datarray),2)]
            dsde  = [float(datarray[i]) for i in range(1,len(datarray),2)]
            
            inter = interpolate.InterpolatedUnivariateSpline(E_lep,dsde)
            PC.E_CC_act = E_1
            if flavor == 0:
                PC.act_dsde_CCe_inter = inter
            elif flavor == 1:
                PC.act_dsde_CCm_inter = inter
            elif flavor == 2:    
                PC.act_dsde_CCt_inter = inter
                
    if flv == 0:
        inter = PC.act_dsde_CCe_inter
        return inter(E_2)[0]
    elif flv == 1:
        inter = PC.act_dsde_CCm_inter
        return inter(E_2)[0]
    elif flv == 2:    
        inter = PC.act_dsde_CCt_inter
        return inter(E_2)[0]
    

def dsde_NC_NusigmaBin_ListOpt(E_1,E_2,param):
    """ Uses the FORTRAN bin to calculate the neutrino/antineutrino
    differential neutral current cross section. The FORTRAN executable
    was generated using the NUSIGMA FORTRAN library. Optimized for
    several calls of the same E_1.
    
    ds/dE in nu + N -> nu + X
    
    Assuming isoscalar target. 
    
    @type   E_1 	:	float
    @param  E_1	    :	neutrino energy [GeV]
    @type   E_2 	:	float
    @param  E_2	    :	neutrino energy [GeV]
    @type   param   :	physics parameter set
    @param  param   :	physics parameter set

    @rtype	        :	float
    @return	        :	neutrino differential NC DIS cross section. New : [eV^-3] Old : [cm^2 GeV^-1]
    """
    if PC.act_dsde_NC_inter == 0 or PC.E_NC_act != E_1 : 
        if param.neutype == "neutrino":
            datarray = subprocess.Popen(global_binpath+'nudsde_sp_list.exe '+str(E_1)+ ' ' +str(1)+' N '+' NC',shell=True,stdout=subprocess.PIPE).communicate()[0]*(param.cm**2/param.GeV)
        elif param.neutype == "antineutrino":
            datarray = subprocess.Popen(global_binpath+'nudsde_sp_list.exe '+str(E_1)+ ' ' +str(2)+' N '+' NC',shell=True,stdout=subprocess.PIPE).communicate()[0]*(param.cm**2/param.GeV)
        else:
            print "Wrong neutrino type"
            quit()
        datarray = datarray.split()
        E_lep = [float(datarray[i]) for i in range(0,len(datarray),2)]
        dsde  = [float(datarray[i]) for i in range(1,len(datarray),2)]
        
        inter = interpolate.InterpolatedUnivariateSpline(E_lep,dsde)
        PC.E_NC_act = E_1
        PC.act_dsde_NC_inter = inter
    else :
        inter = PC.act_dsde_NC_inter
    return inter(E_2)[0]
    
#####################################
## USING nusigma-FORTRAN interface ##
#####################################
    
def dsde_CC_opt(E_1,E_2,param,flv = 0,neutype = None,optimize = 'neutrino'):
    """ Calculates differential neutrino CC DIS cross section. Optimized for
    several calls of the same E_1.
    
    Updated : E_2 optimization mode added.
    
    ds/dE in nu + N -> l + X
    
    Assuming isoscalar target. 
    
    @type   Enu_1 	:	float
    @param  Enu_1	:	neutrino energy [GeV]
    @type   Enu_2 	:	float
    @param  Enu_2	:	lepton energy [GeV]
    @type   param       :	physics parameter set
    @param  param       :	physics parameter set
    @type   flv         :	integer
    @param  flv         :	0 : e, 1 : mu, 2 : tau
    @type   neutype     :	integer
    @param  neutype     :	0 : neutrino, 1 : antineutrino; Overwrites param.neutype.
    @type   neutype     :	string
    @param  neutype     :	neutrino : optimize for E_1, lepton : optimize for E_2
    

    @rtype	        :	float
    @return	        :	neutrino differential CC DIS cross section. New : [eV^-3] Old : [cm^2 GeV^-1]
    """
    if NUDSDE_ERROR :
        quit()
        print "NC:NEU:XSECTIONS:ERROR: Loading NUSIGMA interface : nudsde."
    
    if optimize == "neutrino":
        gen_inter = PC.act_dsde_CCe_n_inter == 0 or PC.act_dsde_CCm_n_inter == 0 or PC.act_dsde_CCt_n_inter == 0 or PC.act_dsde_CCe_a_inter == 0 or PC.act_dsde_CCm_a_inter == 0 or PC.act_dsde_CCt_a_inter == 0 or PC.E_CC_act != E_1
    elif optimize == "lepton":
        gen_inter = PC.act_dsde_CCe_n_inter == 0 or PC.act_dsde_CCm_n_inter == 0 or PC.act_dsde_CCt_n_inter == 0 or PC.act_dsde_CCe_a_inter == 0 or PC.act_dsde_CCm_a_inter == 0 or PC.act_dsde_CCt_a_inter == 0 or PC.E_CC_act != E_2
    
    if gen_inter :
        for flavor in [0,1,2]:
            if optimize == "neutrino":
                E_lep = gt.LogSpaceEnergies(0.1,E_1,200)    # [GeV]
                
                dsde_n = [oxs.dsde(E_1,EE,2*flavor+1,'N','CC')*(param.cm**2/param.GeV) for EE in E_lep]
                dsde_a = [oxs.dsde(E_1,EE,2*flavor+2,'N','CC')*(param.cm**2/param.GeV) for EE in E_lep]
                
                inter_n = interpolate.interp1d(E_lep,dsde_n)
                inter_a = interpolate.interp1d(E_lep,dsde_a)
                
                PC.E_CC_act = E_1
            elif optimize == "lepton":
                E_neu = gt.LogSpaceEnergies(E_2,10000,200)  # [GeV]
                
                dsde_n = [oxs.dsde(EE,E_2,2*flavor+1,'N','CC')*(param.cm**2/param.GeV) for EE in E_neu]
                dsde_a = [oxs.dsde(EE,E_2,2*flavor+2,'N','CC')*(param.cm**2/param.GeV) for EE in E_neu]
                
                inter_n = interpolate.interp1d(E_neu,dsde_n)
                inter_a = interpolate.interp1d(E_neu,dsde_a)
                
                PC.E_CC_act = E_2
            
            if flavor == 0:
                PC.act_dsde_CCe_n_inter = inter_n
                PC.act_dsde_CCe_a_inter = inter_a
            elif flavor == 1:
                PC.act_dsde_CCm_n_inter = inter_n
                PC.act_dsde_CCm_a_inter = inter_a
            elif flavor == 2:    
                PC.act_dsde_CCt_n_inter = inter_n
                PC.act_dsde_CCt_a_inter = inter_a
                
    if neutype == None:
        if param.neutype == "neutrino":
            if flv == 0:
                inter = PC.act_dsde_CCe_n_inter
            elif flv == 1:
                inter = PC.act_dsde_CCm_n_inter
            elif flv == 2:    
                inter = PC.act_dsde_CCt_n_inter
        elif param.neutype == "antineutrino":
            if flv == 0:
                inter = PC.act_dsde_CCe_a_inter
            elif flv == 1:
                inter = PC.act_dsde_CCm_a_inter
            elif flv == 2:    
                inter = PC.act_dsde_CCt_a_inter
    else :
        if neutype == 0:
            if flv == 0:
                inter = PC.act_dsde_CCe_n_inter
            elif flv == 1:
                inter = PC.act_dsde_CCm_n_inter
            elif flv == 2:    
                inter = PC.act_dsde_CCt_n_inter
        elif neutype == 1:
            if flv == 0:
                inter = PC.act_dsde_CCe_a_inter
            elif flv == 1:
                inter = PC.act_dsde_CCm_a_inter
            elif flv == 2:    
                inter = PC.act_dsde_CCt_a_inter
                
    if optimize == "neutrino":
        return inter(E_2)
    elif optimize == "lepton":
        if E_1 == E_2:
            E_1 = E_2+0.00001
        return inter(E_1)
    

def dsde_NC_opt(E_1,E_2,param):
    """ Calculates differential neutrino CC DIS cross section.Optimized for
    several calls of the same E_1.
    
    ds/dE in nu + N -> nu + X
    
    Assuming isoscalar target. 
    
    @type   Enu_1 	:	float
    @param  Enu_1	:	neutrino energy [GeV]
    @type   Enu_2 	:	float
    @param  Enu_2	:	neutrino energy [GeV]
    @type   param       :	physics parameter set
    @param  param       :	physics parameter set

    @rtype	        :	float
    @return	        :	neutrino differential CC DIS cross section. New : [eV^-3] Old : [cm^2 GeV^-1]
    """
    if NUDSDE_ERROR :
        quit()
        print "NC:NEU:XSECTIONS:ERROR: Loading NUSIGMA interface : nudsde."
    
    if PC.act_dsde_NC_n_inter == 0 or PC.act_dsde_NC_a_inter == 0 or PC.E_NC_act != E_1 :
        E_lep = gt.LogSpaceEnergies(0.1,E_1,200)
            
        dsde_n = [oxs.dsde(E_1,EE,1,'N','NC')*(param.cm**2/param.GeV) for EE in E_lep]
        dsde_a = [oxs.dsde(E_1,EE,2,'N','NC')*(param.cm**2/param.GeV) for EE in E_lep]
        
        inter_n = interpolate.interp1d(E_lep,dsde_n)
        inter_a = interpolate.interp1d(E_lep,dsde_a)
        
        PC.E_NC_act = E_1
        PC.act_dsde_NC_n_inter = inter_n
        PC.act_dsde_NC_a_inter = inter_a
        
    if param.neutype == "neutrino":
        inter = PC.act_dsde_NC_n_inter
    elif param.neutype == "antineutrino":
        inter = PC.act_dsde_NC_a_inter
        
    return inter(E_2)
    
#===============================================================================#
# Total Cross Sections                                                          #
#===============================================================================#    
    
#####################################
##       DIS Cross Sections        ##
#####################################
    
def nuDISxsection_CC_NusigmaInt(E,flv,param,datapath = global_datapath,neutype = None):
    """ Returns the total neutrino CC DIS cross section from the dsde
    nusigma cross section.
    
    NOTE :It uses previously generated data from #L{Gen_NuDISxsection_NusigmaInt} using
    the nusigma interface.
    
    sig(nu + N -> lep + X)
    
    Assuming isoscalar target. 
    
    @type   E   	    :	float
    @param  E   	    :	neutrino energy [GeV]
    @type   flv         :	integer
    @param  flv         :	0 : e, 1 : mu, 2 : tau    
    @type   param       :	physics parameter set
    @param  param       :	physics parameter set
    @type   datapath    :	string
    @param  datapath    :	Path to the data
    @type   neutype     :	integer
    @param  neutype     :	0 : neutrino, 1 : antineutrino

    @rtype	        :	float
    @return	        :	neutrino CC DIS cross section. New : [eV^-2] Old : [cm^2]
    """
    if PC.act_sig_CCe_n_inter == 0 or PC.act_sig_CCm_n_inter == 0 or PC.act_sig_CCt_n_inter == 0 or PC.act_sig_CCe_a_inter == 0 or PC.act_sig_CCm_a_inter == 0 or PC.act_sig_CCt_a_inter == 0:
        for flavor in [0,1,2]:
            
            neuneu = param.neutype
            
            filename = "integrated_sig_CC_ineu_"+str(flavor)+"_neutrino.dat"
            file = open(datapath+filename,'r')
            sig_array = []
            gt.hreadfilev4(file,sig_array,param,header_read = False)
            sig_array = sig_array[0]
            param.neutype = neuneu
            
            E_nu  = map(lambda x : x[0],sig_array)
            nxs   = map(lambda x : x[1]*pc.cm**2,sig_array)
            
            filename = "integrated_sig_CC_ineu_"+str(flavor)+"_antineutrino.dat"
            file = open(datapath+filename,'r')
            sig_array = []
            gt.hreadfilev4(file,sig_array,param,header_read = False)
            sig_array = sig_array[0]
            param.neutype = neuneu
            
            E_anu  = map(lambda x : x[0],sig_array)
            axs   = map(lambda x : x[1]*pc.cm**2,sig_array)
                    
            inter_n = interpolate.interp1d(E_nu,nxs)
            inter_a = interpolate.interp1d(E_anu,axs)
        
            if flavor == 0:
                PC.act_sig_CCe_n_inter = inter_n
                PC.act_sig_CCe_a_inter = inter_a
            elif flavor == 1:
                PC.act_sig_CCm_n_inter = inter_n
                PC.act_sig_CCm_a_inter = inter_a
            elif flavor == 2:    
                PC.act_sig_CCt_n_inter = inter_n
                PC.act_sig_CCt_a_inter = inter_a
    if E > 1.0 :
        if neutype == None : 
            if param.neutype == "neutrino":
                if flv == 0:
                    inter = PC.act_sig_CCe_n_inter
                    return inter(E)
                elif flv == 1:
                    inter = PC.act_sig_CCm_n_inter
                    return inter(E)
                elif flv == 2:    
                    inter = PC.act_sig_CCt_n_inter
                    return inter(E)
            elif param.neutype == "antineutrino":
                if flv == 0:
                    inter = PC.act_sig_CCe_a_inter
                    return inter(E)
                elif flv == 1:
                    inter = PC.act_sig_CCm_a_inter
                    return inter(E)
                elif flv == 2:    
                    inter = PC.act_sig_CCt_a_inter
                    return inter(E)
        else :
            if neutype == 0:
                if flv == 0:
                    inter = PC.act_sig_CCe_n_inter
                    return inter(E)
                elif flv == 1:
                    inter = PC.act_sig_CCm_n_inter
                    return inter(E)
                elif flv == 2:    
                    inter = PC.act_sig_CCt_n_inter
                    return inter(E)
            elif neutype == 1:
                if flv == 0:
                    inter = PC.act_sig_CCe_a_inter
                    return inter(E)
                elif flv == 1:
                    inter = PC.act_sig_CCm_a_inter
                    return inter(E)
                elif flv == 2:    
                    inter = PC.act_sig_CCt_a_inter
                    return inter(E)
    else :
        return 0.0
    
def nuDISxsection_NC_NusigmaInt(E,param,datapath = global_datapath,neutype = None):
    """ Returns the total neutrino NC DIS cross section from the dsde
    nusigma cross section.
    
    sig(nu + N -> nu + X)
    
    Assuming isoscalar target. 
    
    NOTE :It uses previously generated data from #L{Gen_NuDISxsection_NusigmaInt} using
    the nusigma interface.    
    
    @type   E   	    :	float
    @param  E   	    :	neutrino energy [GeV]
    @type   flv         :	integer
    @param  flv         :	0 : e, 1 : mu, 2 : tau    
    @type   param       :	physics parameter set
    @param  param       :	physics parameter set
    @type   datapath    :	string
    @param  datapath    :	Path to the data
    @type   neutype     :	integer
    @param  neutype     :	0 : neutrino, 1 : antineutrino

    @rtype	        :	float
    @return	        :	neutrino CC DIS cross section. New : [eV^-2] Old : [cm^2]
    """
    if PC.act_sig_NC_n_inter == 0 or PC.act_sig_NC_a_inter == 0:
        
        neuneu = param.neutype
        
        filename = "integrated_sig_NC_ineu_neutrino.dat"
        file = open(datapath+filename,'r')
        sig_array = []
        gt.hreadfilev4(file,sig_array,param,header_read = False)
        sig_array = sig_array[0]
        param.neutype = neuneu
        
        E_nu  = map(lambda x : x[0],sig_array)
        nxs   = map(lambda x : x[1]*pc.cm**2,sig_array)
        
        filename = "integrated_sig_NC_ineu_antineutrino.dat"
        file = open(datapath+filename,'r')
        sig_array = []
        gt.hreadfilev4(file,sig_array,param,header_read = False)
        sig_array = sig_array[0]
        param.neutype = neuneu
        
        E_anu  = map(lambda x : x[0],sig_array)
        axs   = map(lambda x : x[1]*pc.cm**2,sig_array)
                
        inter_n = interpolate.interp1d(E_nu,nxs)
        inter_a = interpolate.interp1d(E_anu,axs)
        
        PC.act_sig_NC_n_inter = inter_n
        PC.act_sig_NC_a_inter = inter_a
    
    if neutype == None:
        if param.neutype == "neutrino":
            inter = PC.act_sig_NC_n_inter
        elif param.neutype == "antineutrino":
            inter = PC.act_sig_NC_a_inter
    else :
        if neutype == 0:
            inter = PC.act_sig_NC_n_inter
        elif neutype == 1:
            inter = PC.act_sig_NC_a_inter
        
    if E > 1.0:
        return inter(E)
    else : 
        return 0.0

def nuDISxsection_CC_Tbl(Enu,neu):
    """ Returns the total neutrino CC DIS cross section by interpolating
    the tables presented in this reference : 
    
    Ref : Arxiv hep-ph/9807264v1
    
    Also the value for 1GeV has been extrapolated.
    
    @type   Enu   	    :	float
    @param  Enu   	    :	neutrino energy [GeV]
    @type   neu         :	integer
    @param  neu         :	0 : neutrino, 1 : antineutrino

    @rtype	            :	float
    @return	            :	neutrino CC DIS cross section. New : [eV^-2] Old : [cm^2]
    """
    # MANUAL EXTRAPOLATION TO 1 GeV for CC
    m = (0.1932e-36 - 0.7988e-37)/(2.5e1-1.0e1)
    sigCC_1GeV = 0.7988e-37 + (1.0e0-1.0e1)*m

    if Enu < 1.0e0 :
        return 0.0
    else : 
        if(neu == 0):
            E = [1.0e0,1.0e1,2.5e1,6.0e1,1.0e2,2.5e2,6.0e2,1.0e3,2.5e3,6.0e3,1.0e4,2.5e4,6.0e4,1.0e5,2.5e5,6.0e5,1.0e6,2.5e6,6.0e6,1.0e7,2.5e7,6.0e7,1.0e8,2.5e8,6.0e8,1.0e9,2.5e9,6.0e9,1.0e10,2.5e10,6.0e10,1.0e11,2.5e11,6.0e11,1.0e12]
            sigCC = np.array([sigCC_1GeV,0.7988e-37,0.1932e-36,0.4450e-36,0.7221e-36,0.1728e-35,0.3964e-35,0.6399e-35,0.1472e-34,0.3096e-34,0.4617e-34,0.8824e-34,0.1514e-33,0.2022e-33,0.3255e-33,0.4985e-33,0.6342e-33,0.9601e-33,0.1412e-32,0.1749e-32,0.2554e-32,0.3630e-32,0.4436e-32,0.6283e-32,0.8699e-32,0.1049e-31,0.1466e-31,0.2010e-31,0.2379e-31,0.3289e-31,0.4427e-31,0.5357e-31,0.7320e-31,0.9927e-31,0.1179e-30])
            inter=interpolate.interp1d(E,sigCC*pc.cm**2)
        elif(neu == 1):
            E = [1.0e0,1.0e1,2.5e1,6.0e1,1.0e2,2.5e2,6.0e2,1.0e3,2.5e3,6.0e3,1.0e4,2.5e4,6.0e4,1.0e5,2.5e5,6.0e5,1.0e6,2.5e6,6.0e6,1.0e7,2.5e7,6.0e7,1.0e8,2.5e8,6.0e8,1.0e9,2.5e9,6.0e9,1.0e10,2.5e10,6.0e10,1.0e11,2.5e11,6.0e11,1.0e12]
            sigCC = np.array([0.0,0.3936e-37,0.9726e-37,0.2287e-36,0.3747e-36,0.9154e-36,0.2153e-35,0.3542e-35,0.8548e-35,0.1922e-34,0.3008e-34,0.6355e-34,0.1199e-33,0.1683e-33,0.2909e-33,0.4667e-33,0.6051e-33,0.9365e-33,0.1393e-32,0.1734e-32,0.2542e-32,0.3622e-32,0.4430e-32,0.6278e-32,0.8696e-32,0.1050e-31,0.1464e-31,0.2011e-31,0.2406e-31,0.3286e-31,0.4481e-31,0.5335e-31,0.7306e-31,0.9854e-31,0.1165e-30])
            inter=interpolate.interp1d(E,sigCC*pc.cm**2)
        else:
            print "Invalid cross section neutrino type."
            quit()
        return inter(Enu)

def nuDISxsection_NC_Tbl(Enu,neu):
    """ Returns the total neutrino NC DIS cross section by interpolating
    the tables presented in this reference : 
    
    Ref : Arxiv hep-ph/9807264v1
    
    Also the value for 1GeV has been extrapolated.
    
    @type   Enu   	    :	float
    @param  Enu   	    :	neutrino energy [GeV]
    @type   neu         :	integer
    @param  neu         :	0 : neutrino, 1 : antineutrino

    @rtype	            :	float
    @return	            :	neutrino NC DIS cross section. New : [eV^-2] Old : [cm^2]
    """
    if Enu < 1.0e1 :
        return 0.0
    else :
        if(neu == 0):
            E = [1.0e1,2.5e1,6.0e1,1.0e2,2.5e2,6.0e2,1.0e3,2.5e3,6.0e3,1.0e4,2.5e4,6.0e4,1.0e5,2.5e5,6.0e5,1.0e6,2.5e6,6.0e6,1.0e7,2.5e7,6.0e7,1.0e8,2.5e8,6.0e8,1.0e9,2.5e9,6.0e9,1.0e10,2.5e10,6.0e10,1.0e11,2.5e11,6.0e11,1.0e12]
            sigNC = np.array([0.2492e-37,0.6033e-37,0.1391e-36,0.2261e-36,0.5430e-36,0.1255e-36,0.2039e-35,0.4781e-35,0.1035e-34,0.1575e-34,0.3139e-34,0.5615e-34,0.7667e-34,0.1280e-33,0.2017e-33,0.2600e-33,0.4018e-33,0.6001e-33,0.7482e-33,0.1104e-32,0.1581e-32,0.1939e-32,0.2763e-32,0.3837e-32,0.4641e-32,0.6490e-32,0.8931e-32,0.1066e-31,0.1465e-31,0.1995e-31,0.2377e-31,0.3247e-31,0.4377e-31,0.5196e-31])
            inter=interpolate.interp1d(E,sigNC*pc.cm**2 )
        elif(neu == 1):
            E = [1.0e1,2.5e1,6.0e1,1.0e2,2.5e2,6.0e2,1.0e3,2.5e3,6.0e3,1.0e4,2.5e4,6.0e4,1.0e5,2.5e5,6.0e5,1.0e6,2.5e6,6.0e6,1.0e7,2.5e7,6.0e7,1.0e8,2.5e8,6.0e8,1.0e9,2.5e9,6.0e9,1.0e10,2.5e10,6.0e10,1.0e11,2.5e11,6.0e11,1.0e12]
            sigNC = np.array([0.1381e-37,0.3403e-37,0.7982e-37,0.1307e-36,0.3193e-36,0.7531e-36,0.1243e-35,0.3026e-35,0.6896e-35,0.1091e-34,0.2358e-34,0.4570e-34,0.6515e-34,0.1158e-33,0.1901e-33,0.2493e-33,0.3929e-33,0.5930e-33,0.7423e-33,0.1100e-32,0.1578e-32,0.1937e-32,0.2762e-32,0.3836e-32,0.4641e-32,0.6489e-32,0.8931e-32,0.1066e-31,0.1465e-31,0.1995e-31,0.2377e-31,0.3247e-31,0.4377e-31,0.5195e-31])
            inter=interpolate.interp1d(E,sigNC*pc.cm**2 )
        else:
            print "Invalid cross section neutrino type."
            quit()
        return inter(Enu)

def nuDISxsection_NCANDCC_Tbl(Enu,neu):
    """ Returns the total neutrino NC+CC DIS cross section by interpolating
    the tables presented in this reference : 
    
    Ref : Arxiv hep-ph/9807264v1
    
    Also the value for 1GeV has been extrapolated.
    
    @type   Enu   	    :	float
    @param  Enu   	    :	neutrino energy [GeV]
    @type   neu         :	integer
    @param  neu         :	0 : neutrino, 1 : antineutrino

    @rtype	            :	float
    @return	            :	neutrino NC+CC DIS cross section. New : [eV^-2] Old : [cm^2]
    """
    if Enu < 1.0e1 :
        return 0.0
    else : 
        if(neu == 0):
            E = [1.0e1,2.5e1,6.0e1,1.0e2,2.5e2,6.0e2,1.0e3,2.5e3,6.0e3,1.0e4,2.5e4,6.0e4,1.0e5,2.5e5,6.0e5,1.0e6,2.5e6,6.0e6,1.0e7,2.5e7,6.0e7,1.0e8,2.5e8,6.0e8,1.0e9,2.5e9,6.0e9,1.0e10,2.5e10,6.0e10,1.0e11,2.5e11,6.0e11,1.0e12]
            sigtot = np.array([0.1048e-36,0.2535e-36,0.5841e-36,0.9482e-36,0.2271e-35,0.5219e-35,0.8438e-35,0.1950e-34,0.4131e-34,0.6192e-34,0.1196e-33,0.2076e-33,0.2789e-33,0.4535e-33,0.7002e-33,0.8942e-33,0.1362e-32,0.2012e-32,0.2497e-32,0.3658e-32,0.5211e-32,0.6375e-32,0.9046e-32,0.1254e-31,0.1513e-31,0.2115e-31,0.2903e-13,0.3445e-31,0.4754e-31,0.6422e-31,0.7734e-31,0.1057e-30,0.1430e-30,0.1699e-30])
            inter=interpolate.interp1d(E,sigtot*pc.cm**2 )
        elif(neu == 1):
            E = [1.0e1,2.5e1,6.0e1,1.0e2,2.5e2,6.0e2,1.0e3,2.5e3,6.0e3,1.0e4,2.5e4,6.0e4,1.0e5,2.5e5,6.0e5,1.0e6,2.5e6,6.0e6,1.0e7,2.5e7,6.0e7,1.0e8,2.5e8,6.0e8,1.0e9,2.5e9,6.0e9,1.0e10,2.5e10,6.0e10,1.0e11,2.5e11,6.0e11,1.0e12]
            sigtot = np.array([0.5317e-37,0.1313e-36,0.3085e-36,0.5054e-36,0.1235e-35,0.2906e-35,0.4785e-35,0.1157e-34,0.2612e-34,0.4099e-34,0.8713e-34,0.1656e-33,0.2334e-33,0.4067e-33,0.6568e-33,0.8544e-33,0.1329e-32,0.1986e-32,0.2476e-32,0.3642e-32,0.5200e-32,0.6367e-32,0.9040e-32,0.1253e-31,0.1514e-31,0.2113e-31,0.2904e-31,0.3472e-31,0.4751e-31,0.6476e-31,0.7712e-31,0.1055e-30,0.1423e-30,0.1685e-30])
            inter=interpolate.interp1d(E,sigtot*pc.cm**2)
        else:
            print "Invalid cross section neutrino type."
            quit()
        return inter(Enu)
		
#####################################
##    Low energy cross Sections    ##
#####################################		

def nuDeuteriumxsection_CC_Tbl(Enu,neu,return_interpolator = False):
    """ Returns the electron-neutrino/antineutrino-Deuterium CC DIS cross section by interpolating
    the tables presented in this reference : 
    
    Ref : Int. J. Mod. Phys. E 3, 101 (1994)
    
    @type   Enu   	    :	float
    @param  Enu   	    :	neutrino energy [GeV]
    @type   neu         :	integer
    @param  neu         :	0 : neutrino, 1 : antineutrino
    @type   return_interpolator  :	boolean
    @param  return_interpolator  :	if True/False : will return interolator/value

    @rtype	            :	float
    @return	            :	neutrino CC DIS cross section. New : [eV^-2] Old : [cm^2]
    """
    
    if neu == 0:
        print "Missing nu-Deuterium cross section."
    elif neu == 1:    
        # NOTE : here the cross section was given in units of 10^{-42}
        E1 = np.array(np.append(np.append(np.append(np.append(np.arange(2.0,12.2,0.2),np.arange(12.5,25.5,0.5)),np.arange(26.0,56.0,1.0)),np.arange(60.0,105.0,5.0)),np.arange(110.0,180.0,10.0)))
        sigCC1 = np.array([1.e-30,1.e-30,1.e-30,1.e-30,1.e-30,1.e-30,1.e-30,1.e-30,1.e-30,1.e-30,1.e-30,1.065e-3,4.397e-3,9.832e-3,1.748e-2,2.747e-2,3.993e-2,5.498e-2,7.270e-2,9.319e-2,1.165e-1,1.427e-1,1.719e-1,2.040e-1,2.392e-1,2.774e-1,3.187e-1,3.631e-1,4.106e-1,4.612e-1,5.150e-1,5.719e-1,6.320e-1,6.967e-1,7.635e-1,8.336e-1,9.068e-1,9.833e-1,1.063,1.146,1.232,1.322,1.415,1.511,1.610,1.712,1.818,1.927,2.039,2.155,2.273,2.584,2.914,3.265,3.635,4.025,4.435,4.864,5.313,5.780,6.267,6.773,7.297,7.841,8.402,8.983,9.581,1.020e1,1.083e1,1.149e1,1.216e1,1.285e1,1.355e1,1.427e1,1.501e1,1.577e1,1.655e1,1.815e1,1.981e1,2.155e1,2.335e1,2.521e1,2.713e1,2.912e1,3.117e1,3.328e1,3.545e1,3.768e1,3.996e1,4.231e1,4.470e1,4.716e1,4.966e1,5.222e1,5.484e1,5.750e1,6.021e1,6.298e1,6.579e1,6.865e1,7.156e1,7.451e1,7.751e1,8.055e1,8.364e1,8.676e1,8.993e1,1.064e2,1.237e2,1.419e2,1.607e2,1.801e2,2.001e2,2.205e2,2.412e2,2.623e2,3.049e2,3.481e2,3.914e2,4.346e2,4.776e2,5.201e2,5.623e2])           
        inter=interpolate.interp1d(E1,sigCC1*1.0e-42*pc.cm**2)
        if return_interpolador :
            return inter
        else :
            return inter(Enu)
    else : 
        print "Invalid cross section neutrino type."

def nuDeuteriumxsection_NC_Tbl(Enu,neu,return_interpolator = False):
    """ Returns the electron-neutrino/antineutrino-Deuterium NC DIS cross section by interpolating
    the tables presented in this reference : 
    
    Ref : Int. J. Mod. Phys. E 3, 101 (1994)
    
    @type   Enu   	    :	float
    @param  Enu   	    :	neutrino energy [GeV]
    @type   neu         :	integer
    @param  neu         :	0 : neutrino, 1 : antineutrino

    @rtype	            :	float
    @return	            :	neutrino NC DIS cross section. New : [eV^-2] Old : [cm^2]
    """
    if neu == 0:
        print "Missing nu-Deuterium cross section."
    elif neu == 1:    
        # NOTE : here the cross section was given in units of 10^{-42}
        E1 = np.array(np.append(np.append(np.append(np.append(np.arange(2.0,12.2,0.2),np.arange(12.5,25.5,0.5)),np.arange(26.0,56.0,1.0)),np.arange(60.0,105.0,5.0)),np.arange(110.0,180.0,10.0)))
        sigNC1 = np.array([1.e-30,1.e-30,4.362e-5,4.253e-4,1.451e-3,3.334e-3,6.236e-3,1.028e-2,1.557e-2,2.219e-2,3.021e-2,3.967e-2,5.064e-2,6.314e-2,7.722e-2,9.290e-2,1.102e-1,1.292e-1,1.498e-1,1.721e-1,1.961e-1,2.218e-1,2.491e-1,2.782e-1,3.102e-1,3.430e-1,3.776e-1,4.140e-1,4.522e-1,4.921e-1,5.339e-1,5.775e-1,6.228e-1,6.700e-1,7.189e-1,7.697e-1,8.223e-1,8.767e-1,9.328e-1,9.908e-1,1.051,1.112,1.176,1.241,1.308,1.377,1.447,1.520,1.594,1.670,1.748,1.950,2.164,2.389,2.625,2.872,3.130,3.399,3.679,3.969,4.271,4.584,4.907,5.241,5.585,5.940,6.306,6.682,7.069,7.466,7.873,8.291,8.719,9.157,9.606,1.006e1,1.053e1,1.150e1,1.251e1,1.355e1,1.464e1,1.576e1,1.692e1,1.812e1,1.936e1,2.063e1,2.194e1,2.329e1,2.467e1,2.609e1,2.754e1,2.902e1,3.055e1,3.210e1,3.368e1,3.530e1,3.695e1,3.864e1,4.035e1,4.210e1,4.387e1,4.568e1,4.751e1,4.937e1,5.126e1,5.318e1,5.513e1,6.525e1,7.599e1,8.728e1,9.906e1,1.113e2,1.239e2,1.368e2,1.500e2,1.634e2,1.908e2,2.186e2,2.465e2,2.744e2,3.020e2,3.292e2,3.559e2])
        inter=interpolate.interp1d(E1,sigNC1*1.0e-42*pc.cm**2)
        if return_interpolador :
            return inter
        else :
            return inter(Enu)
    else : 
        print "Invalid cross section neutrino type."  

def nuDeuteriumxsection_NCANDCC_Tbl(Enu,neu,return_interpolator = False):
    """ Returns the electron-neutrino/antineutrino-Deuterium total (NC+CC) DIS cross section by interpolating
    the tables presented in this reference : 
    
    Ref : Int. J. Mod. Phys. E 3, 101 (1994)
    
    @type   Enu   	    :	float
    @param  Enu   	    :	neutrino energy [GeV]
    @type   neu         :	integer
    @param  neu         :	0 : neutrino, 1 : antineutrino

    @rtype	            :	float
    @return	            :	neutrino NC+CC DIS cross section. New : [eV^-2] Old : [cm^2]
    """
    if neu == 0:
        print "Missing nu-Deuterium cross section."
    elif neu == 1:    
        # NOTE : here the cross section was given in units of 10^{-42}
        E1 = np.array(np.append(np.append(np.append(np.append(np.arange(2.0,12.2,0.2),np.arange(12.5,25.5,0.5)),np.arange(26.0,56.0,1.0)),np.arange(60.0,105.0,5.0)),np.arange(110.0,180.0,10.0)))
        sigCC1 = np.array([1.e-30,1.e-30,1.e-30,1.e-30,1.e-30,1.e-30,1.e-30,1.e-30,1.e-30,1.e-30,1.e-30,1.065e-3,4.397e-3,9.832e-3,1.748e-2,2.747e-2,3.993e-2,5.498e-2,7.270e-2,9.319e-2,1.165e-1,1.427e-1,1.719e-1,2.040e-1,2.392e-1,2.774e-1,3.187e-1,3.631e-1,4.106e-1,4.612e-1,5.150e-1,5.719e-1,6.320e-1,6.967e-1,7.635e-1,8.336e-1,9.068e-1,9.833e-1,1.063,1.146,1.232,1.322,1.415,1.511,1.610,1.712,1.818,1.927,2.039,2.155,2.273,2.584,2.914,3.265,3.635,4.025,4.435,4.864,5.313,5.780,6.267,6.773,7.297,7.841,8.402,8.983,9.581,1.020e1,1.083e1,1.149e1,1.216e1,1.285e1,1.355e1,1.427e1,1.501e1,1.577e1,1.655e1,1.815e1,1.981e1,2.155e1,2.335e1,2.521e1,2.713e1,2.912e1,3.117e1,3.328e1,3.545e1,3.768e1,3.996e1,4.231e1,4.470e1,4.716e1,4.966e1,5.222e1,5.484e1,5.750e1,6.021e1,6.298e1,6.579e1,6.865e1,7.156e1,7.451e1,7.751e1,8.055e1,8.364e1,8.676e1,8.993e1,1.064e2,1.237e2,1.419e2,1.607e2,1.801e2,2.001e2,2.205e2,2.412e2,2.623e2,3.049e2,3.481e2,3.914e2,4.346e2,4.776e2,5.201e2,5.623e2])                
        sigNC1 = np.array([1.e-30,1.e-30,4.362e-5,4.253e-4,1.451e-3,3.334e-3,6.236e-3,1.028e-2,1.557e-2,2.219e-2,3.021e-2,3.967e-2,5.064e-2,6.314e-2,7.722e-2,9.290e-2,1.102e-1,1.292e-1,1.498e-1,1.721e-1,1.961e-1,2.218e-1,2.491e-1,2.782e-1,3.102e-1,3.430e-1,3.776e-1,4.140e-1,4.522e-1,4.921e-1,5.339e-1,5.775e-1,6.228e-1,6.700e-1,7.189e-1,7.697e-1,8.223e-1,8.767e-1,9.328e-1,9.908e-1,1.051,1.112,1.176,1.241,1.308,1.377,1.447,1.520,1.594,1.670,1.748,1.950,2.164,2.389,2.625,2.872,3.130,3.399,3.679,3.969,4.271,4.584,4.907,5.241,5.585,5.940,6.306,6.682,7.069,7.466,7.873,8.291,8.719,9.157,9.606,1.006e1,1.053e1,1.150e1,1.251e1,1.355e1,1.464e1,1.576e1,1.692e1,1.812e1,1.936e1,2.063e1,2.194e1,2.329e1,2.467e1,2.609e1,2.754e1,2.902e1,3.055e1,3.210e1,3.368e1,3.530e1,3.695e1,3.864e1,4.035e1,4.210e1,4.387e1,4.568e1,4.751e1,4.937e1,5.126e1,5.318e1,5.513e1,6.525e1,7.599e1,8.728e1,9.906e1,1.113e2,1.239e2,1.368e2,1.500e2,1.634e2,1.908e2,2.186e2,2.465e2,2.744e2,3.020e2,3.292e2,3.559e2])
        inter=interpolate.interp1d(E1,(sigNC1+sigCC1)*1.0e-42*pc.cm**2)
        if return_interpolador :
            return inter
        else :
            return inter(Enu)
    else : 
        print "Invalid cross section neutrino type."
    
def nuFexsection_CC_Tbl(Enu,neu,return_interpolator = False):
    """ Returns electron-the neutrino/antineutrino-Fe(56) CC DIS cross section by interpolating
    the tables presented in this reference : 
    
    Ref : Phys. Rev. C 63, 025802 (2001)
    
    @type   Enu   	    :	float
    @param  Enu   	    :	neutrino energy [GeV]
    @type   neu         :	integer
    @param  neu         :	0 : neutrino, 1 : antineutrino
    @type   return_interpolator  :	boolean
    @param  return_interpolator  :	if True/False : will return interolator/value

    @rtype	            :	float
    @return	            :	neutrino/antineutrino-Fe(56) CC DIS cross section. New : [eV^-2] Old : [cm^2]
    """
    if neu == 0:
        # NOTE : here the cross section was given in units of 10^{-42}
        E2 = np.arange(10.0,155.0,5.0)
        sigCC2 = np.array([6.61e-1,6.45,2.93e1,7.33e1,1.40e2,2.36e2,3.71e2,5.55e2,7.98e2,1.10e3,1.48e3,1.92e3,2.42e3,2.99e3,3.60e3,4.27e3,4.98e3,5.73e3,6.52e3,7.36e3,8.24e3,9.16e3,1.01e4,1.11e4,1.21e4,1.32e4,1.42e4,1.53e4,1.64e4])
        inter=interpolate.interp1d(E2,sigCC2*1.0e-42*pc.cm**2)
        if return_interpolador :
            return inter
        else :
            return inter(Enu)
    elif neu == 1:    
        print "Missing antinu-Pb cross section."
    else : 
        print "Invalid cross section neutrino type."

def nuFexsection_NC_Tbl(Enu,neu,return_interpolator = False):
    """ Returns the electron-neutrino/antineutrino-Fe(56) NC cross section by interpolating
    the tables presented in this reference : 
    
    Ref : Phys. Rev. C 63, 025802 (2001)
    
    @type   Enu   	    :	float
    @param  Enu   	    :	neutrino energy [GeV]
    @type   neu         :	integer
    @param  neu         :	0 : neutrino, 1 : antineutrino

    @rtype	            :	float
    @return	            :	neutrino/antineutrino-Fe(56) NC DIS cross section. New : [eV^-2] Old : [cm^2]
    """
    if neu == 0:
        # NOTE : here the cross section was given in units of 10^{-42}
        E2 = np.arange(10.0,155.0,5.0)
        sigNC2 = np.array([1.91e-1,2.19,6.90,1.51e1,2.85e1,4.89e1,7.86e1,1.19e2,1.72e2,2.39e2,3.20e2,4.15e2,5.25e2,6.50e2,7.89e2,9.42e2,1.11e3,1.29e3,1.49e3,1.70e3,1.92e3,2.16e3,2.41e3,2.66e3,2.92e3,3.19e3,3.46e3,3.74e3,4.01e3])
        inter=interpolate.interp1d(E2,sigNC2*1.0e-42*pc.cm**2)
        if return_interpolador :
            return inter
        else :
            return inter(Enu)
    elif neu == 1:    
        print "Missing antinu-Pb cross section."
    else : 
        print "Invalid cross section neutrino type."

def nuFexsection_NCANDCC_Tbl(Enu,neu,return_interpolator = False):
    """ Returns the electron-neutrino/antineutrino-Fe(56) total cross section by interpolating
    the tables presented in this reference : 
    
    Ref : Phys. Rev. C 63, 025802 (2001)
    
    @type   Enu   	    :	float
    @param  Enu   	    :	neutrino energy [GeV]
    @type   neu         :	integer
    @param  neu         :	0 : neutrino, 1 : antineutrino

    @rtype	            :	float
    @return	            :	neutrino/antineutrino-Fe(56) NC+CC DIS cross section. New : [eV^-2] Old : [cm^2]
    """
    if neu == 0:
        # NOTE : here the cross section was given in units of 10^{-42}
        E2 = np.arange(10.0,155.0,5.0)
        sigCC2 = np.array([6.61e-1,6.45,2.93e1,7.33e1,1.40e2,2.36e2,3.71e2,5.55e2,7.98e2,1.10e3,1.48e3,1.92e3,2.42e3,2.99e3,3.60e3,4.27e3,4.98e3,5.73e3,6.52e3,7.36e3,8.24e3,9.16e3,1.01e4,1.11e4,1.21e4,1.32e4,1.42e4,1.53e4,1.64e4])
        sigNC2 = np.array([1.91e-1,2.19,6.90,1.51e1,2.85e1,4.89e1,7.86e1,1.19e2,1.72e2,2.39e2,3.20e2,4.15e2,5.25e2,6.50e2,7.89e2,9.42e2,1.11e3,1.29e3,1.49e3,1.70e3,1.92e3,2.16e3,2.41e3,2.66e3,2.92e3,3.19e3,3.46e3,3.74e3,4.01e3])
        inter=interpolate.interp1d(E2,(sigCC2+sigNC2)*1.0e-42*pc.cm**2)
        if return_interpolador :
            return inter
        else :
            return inter(Enu)
    elif neu == 1:    
        print "Missing antinu-Pb cross section."
    else : 
        print "Invalid cross section neutrino type."

def nuPbxsection_CC_Tbl(Enu,neu,return_interpolator = False):
    """ Returns the electron-neutrino/antineutrino-Pb(208) CC DIS cross section by interpolating
    the tables presented in this reference : 
    
    Ref : Phys. Rev. C 63, 025802 (2001)
    
    @type   Enu   	    :	float
    @param  Enu   	    :	neutrino energy [GeV]
    @type   neu         :	integer
    @param  neu         :	0 : neutrino, 1 : antineutrino
    @type   return_interpolator  :	boolean
    @param  return_interpolator  :	if True/False : will return interolator/value

    @rtype	            :	float
    @return	            :	neutrino/antineutrino-Pb(208) CC DIS cross section. New : [eV^-2] Old : [cm^2]
    """
    if neu == 0:
        # NOTE : here the cross section was given in units of 10^{-42}
        E3 = np.arange(10.0,155.0,5.0)
        sigCC3 = np.array([9.34,1.41e2,4.85e2,1.32e3,2.48e3,3.99e3,5.72e3,7.63e3,9.69e3,1.20e4,1.45e4,1.73e4,2.02e4,2.31e4,2.62e4,2.93e4,3.26e4,3.60e4,3.96e4,4.33e4,4.71e4,5.10e4,5.50e4,5.90e4,6.31e4,6.71e4,7.12e4,7.52e4,7.91e4])
        inter=interpolate.interp1d(E3,sigCC3*1.0e-42*pc.cm**2)
        if return_interpolador :
            return inter
        else :
            return inter(Enu)
    elif neu == 1:    
        print "Missing antinu-Pb cross section."
    else : 
        print "Invalid cross section neutrino type."

def nuPbxsection_NC_Tbl(Enu,neu,return_interpolator = False):
    """ Returns the electron-neutrino/antineutrino-Pb(208) NC cross section by interpolating
    the tables presented in this reference : 
    
    Ref : Phys. Rev. C 63, 025802 (2001)
    
    @type   Enu   	    :	float
    @param  Enu   	    :	neutrino energy [GeV]
    @type   neu         :	integer
    @param  neu         :	0 : neutrino, 1 : antineutrino

    @rtype	            :	float
    @return	            :	neutrino/antineutrino-Pb(208) NC DIS cross section. New : [eV^-2] Old : [cm^2]
    """
    if neu == 0:
        # NOTE : here the cross section was given in units of 10^{-42}
        E3 = np.arange(10.0,155.0,5.0)
        sigNC3 = np.array([7.14e-1,7.98,2.54e1,5.84e1,1.14e2,1.99e2,3.17e2,4.72e2,6.65e2,8.96e2,1.17e3,1.48e3,1.83e3,2.22e3,2.65e3,3.11e3,3.61e3,4.13e3,4.69e3,5.26e3,5.86e3,6.47e3,7.09e3,7.73e3,8.37e3,9.01e3,9.66e3,1.03e4,1.09e4])
        inter=interpolate.interp1d(E3,sigNC3*1.0e-42*pc.cm**2)
        if return_interpolador :
            return inter
        else :
            return inter(Enu)
    elif neu == 1:    
        print "Missing antinu-Pb cross section."
    else : 
        print "Invalid cross section neutrino type."

def nuPbxsection_NCANDCC_Tbl(Enu,neu,return_interpolator = False):
    """ Returns the electron-neutrino/antineutrino-Pb(208) total cross section by interpolating
    the tables presented in this reference : 
    
    Ref : Phys. Rev. C 63, 025802 (2001)
    
    @type   Enu   	    :	float
    @param  Enu   	    :	neutrino energy [GeV]
    @type   neu         :	integer
    @param  neu         :	0 : neutrino, 1 : antineutrino

    @rtype	            :	float
    @return	            :	neutrino/antineutrino-Pb(208) NC+CC DIS cross section. [cm^2]
    """
    if neu == 0:
        # NOTE : here the cross section was given in units of 10^{-42}
        E3 = np.arange(10.0,155.0,5.0)
        sigCC3 = np.array([9.34,1.41e2,4.85e2,1.32e3,2.48e3,3.99e3,5.72e3,7.63e3,9.69e3,1.20e4,1.45e4,1.73e4,2.02e4,2.31e4,2.62e4,2.93e4,3.26e4,3.60e4,3.96e4,4.33e4,4.71e4,5.10e4,5.50e4,5.90e4,6.31e4,6.71e4,7.12e4,7.52e4,7.91e4])
        sigNC3 = np.array([7.14e-1,7.98,2.54e1,5.84e1,1.14e2,1.99e2,3.17e2,4.72e2,6.65e2,8.96e2,1.17e3,1.48e3,1.83e3,2.22e3,2.65e3,3.11e3,3.61e3,4.13e3,4.69e3,5.26e3,5.86e3,6.47e3,7.09e3,7.73e3,8.37e3,9.01e3,9.66e3,1.03e4,1.09e4])
        inter=interpolate.interp1d(E3,(sigCC3+sigNC3)*1.0e-42*pc.cm**2)
        if return_interpolador :
            return inter
        else :
            return inter(Enu)
    elif neu == 1:    
        print "Missing antinu-Pb cross section."
    else : 
        print "Invalid cross section neutrino type."
        
def numuPbxsection_CC_Tbl(Enu,neu,return_interpolator = False):
    """ Returns the muon-neutrino/antineutrino-Pb(208) CC DIS cross section by interpolating
    the tables presented in this reference : 
    
    Ref : MISSING REFERENCE.
    
    @type   Enu   	    :	float
    @param  Enu   	    :	neutrino energy [GeV]
    @type   neu         :	integer
    @param  neu         :	0 : neutrino, 1 : antineutrino
    @type   return_interpolator  :	boolean
    @param  return_interpolator  :	if True/False : will return interolator/value

    @rtype	            :	float
    @return	            :	neutrino/antineutrino-Pb(208) CC DIS cross section. New : [eV^-2] Old : [cm^2]
    """
    if neu == 0:
        print "Missing nu_mu-Pb cross section."
    elif neu == 1:    
        # NOTE : here the cross section was given in units of 10^{-40} cm^2
        Enu_mu_MeV = [100.0,105.0,110.303467576,115.450609996,118.779570598,121.804189077,123.920591193,125.733963007,127.547772093,129.964580874,131.1732039,132.985263894,134.493419039,136.608946609,138.117101753,140.528663256,142.64331628,144.752722025,148.372469282,149.878000787,153.49731077,156.81227863,160.132056496,167.383357383,171.61528707,185.230224321,192.796361887,200.667278849,207.329135511,217.020420657,229.136822773,246.40736368,259.135948227,265.196554287,284.292273383,293.688399143,274.895273077]
        Enu_mu_GeV = map(lambda x : x*1.0e-3, Enu_mu_MeV)
        Sigmupb = np.array([0.0,0.0,0.0,14.43001443,28.86002886,47.619047619,63.4920634921,79.3650793651,93.7950937951,118.326118326,129.87012987,150.072150072,173.16017316,215.007215007,191.919191919,256.854256854,317.46031746,278.499278499,372.294372294,404.04040404,460.317460317,520.923520924,565.656565657,636.363636364,670.995670996,741.702741703,773.448773449,799.422799423,815.295815296,834.054834055,849.927849928,857.142857143,852.813852814,852.813852814,836.940836941,829.725829726,847.041847042])
        inter=interpolate.interp1d(Enu_mu_GeV,Sigmupb*1.0e-40*pc.cm**2)
        if return_interpolador :
            return inter
        else :
            if Enu < 110.0e-3 :#[GeV]
                return 0.0
            else :
                return inter(Enu)
    else : 
        print "Invalid cross section neutrino type."  
        
def nueCxsection_CC_Tbl(Enu,neu,return_interpolator = False):
    """ Returns the muon-neutrino/antineutrino-C(XXX) CC DIS cross section by interpolating
    the tables presented in this reference : 
    
    Ref : arXiv : 1005.2134
    
    @type   Enu   	    :	float
    @param  Enu   	    :	neutrino energy [GeV]
    @type   neu         :	integer
    @param  neu         :	0 : neutrino, 1 : antineutrino
    @type   return_interpolator  :	boolean
    @param  return_interpolator  :	if True/False : will return interolator/value

    @rtype	            :	float
    @return	            :	neutrino/antineutrino-C(XXX) CC DIS cross section. New : [eV^-2] Old : [cm^2]
    """
    if neu == 0:
        print "Missing nu_mu-Pb cross section."
    elif neu == 1:    
        # NOTE : here the cross section was given in units of 10^{-39} cm^2
        Enu_mu_MeV = np.arange(0.0,600.0+10.0,10.0)
        Enu_mu_GeV = map(lambda x : x*1.0e-3, Enu_mu_MeV)
        SigC = np.array([0.0,0.0,0.01,0.02,0.04,0.08,0.12,0.16,0.22,0.28,0.34,0.42,0.49,0.58,0.67,0.76,0.85,0.95,1.06,1.16,1.27,1.39,1.5,1.62,1.74,1.86,1.99,2.12,2.25,2.38,2.52,2.66,2.8,2.94,3.09,3.23,3.38,3.53,3.69,3.84,3.99,4.15,4.31,4.47,4.62,4.79,4.95,5.11,5.27,5.43,5.6,5.76,5.93,6.09,6.26,6.43,6.59,6.76,6.92,7.09,7.26])
        inter=interpolate.interp1d(Enu_mu_GeV,SigC*1.0e-39*pc.cm**2)
        if return_interpolator :
            return inter
        else :
            return inter(Enu)
    else : 
        print "Invalid cross section neutrino type."  
        
def numuCxsection_CC_Tbl(Enu,neu,return_interpolator = False):
    """ Returns the muon-neutrino/antineutrino-C(XXX) CC DIS cross section by interpolating
    the tables presented in this reference : 
    
    Ref : arXiv : 1005.2134
    
    @type   Enu   	    :	float
    @param  Enu   	    :	neutrino energy [GeV]
    @type   neu         :	integer
    @param  neu         :	0 : neutrino, 1 : antineutrino
    @type   return_interpolator  :	boolean
    @param  return_interpolator  :	if True/False : will return interolator/value

    @rtype	            :	float
    @return	            :	neutrino/antineutrino-C(XXX) CC DIS cross section. New : [eV^-2] Old : [cm^2]
    """
    if neu == 0:
        print "Missing nu_mu-Pb cross section."
    elif neu == 1:    
        # NOTE : here the cross section was given in units of 10^{-39} cm^2
        Enu_mu_GeV = [0.2994232308,0.3524765761,0.4025377843,0.4511064551,0.4996727718,0.5445171618,0.5930952493,0.6356984792,0.6805522859,0.725413155,0.7665332643,0.8121404021,0.8577475399,0.8988747116,0.9437261641,0.9803663073,1.0177550732,1.0543999247,1.1007603936]
        SigC = [0.0340694006,0.0477917981,0.0600946372,0.0723974763,0.0842271293,0.0979495268,0.1121451104,0.1253943218,0.1410094637,0.158044164,0.1731861199,0.1902208202,0.2072555205,0.2238170347,0.2389589905,0.2536277603,0.2687697161,0.284384858,0.3028391167]
        inter=interpolate.interp1d(Enu_mu_GeV,SigC*1.0e-39*pc.cm**2)
        if return_interpolator :
            return inter
        else :
            if Enu < 300.0e-3 :#[GeV]
                return 0.0
            else :
                return inter(Enu)
    else : 
        print "Invalid cross section neutrino type."
		
#===============================================================================#
# Analytic differential cross section calculation - Experimental                #
#===============================================================================# 		
                
def difsigDISNC(Enu_1,Enu_2,track,body,param):
    """ Calculates differential neutrino NC DIS cross section.
    
    d/dE_2_nu sig(\nu_1 N -> \nu_2 X)
    
    Ref.: Appendix B.  arXiv: hep-ph/0506298
    
    @type   Enu_1 	:	float
    @param  Enu_1	:	neutrino energy [GeV]
    @type   Enu_2 	:	float
    @param  Enu_2	:	neutrino energy [GeV]
    @type   track 	:	track object
    @param  track	:	trayectory
    @type   body 	:	body
    @param  body	:	body [Sun,Earth,Vacuum]
    @type   param   :	physics parameter set
    @param  param   :	physics parameter set

    @rtype	        :	float
    @return	        :	neutrino differential NC DIS cross section. New : [eV^-3] Old : [cm^2 eV^-1]
    """
    
    # Z-coupling fo quarks
    
    g_Lu = 0.5 -(2.0/3.0)*param.sw_sq
    g_Ru = -(2.0/3.0)*param.sw_sq
    g_Ld = -0.5 +(2.0/3.0)**param.sw_sq
    g_Rd = (1.0/3.0)*param.sw_sq
    
    # neutron and proton number density
    
    nucleon_density  = body.density(track)*param.Na # [nucleon/cm^3]
    ye               = body.ye(track)
    Np,Nn            = nucleon_density*ye,nucleon_density*(1-ye)
    
    # quark transfer momentum ratio
    
    p_u     = (0.25*Np + 0.15*Nn)/(Np+Nn)
    p_d     = (0.25*Nn + 0.15*Np)/(Np+Nn)
    
    p_u_b   = (0.03*Np + 0.06*Nn)/(Np+Nn)
    p_d_b   = (0.03*Nn + 0.06*Np)/(Np+Nn)
    
    m_n     = param.neutron_mass*param.GeV
    
    # calculatin xsections
    
    if param.neutype == "neutrino":
        dsigdE_u = ((2.0*param.GF**2*m_n)/np.pi)*(p_u*(g_Lu**2+g_Ru**2*(Enu_2/Enu_1)**2)+p_u_b*(g_Ru**2+g_Lu**2*(Enu_2/Enu_1)**2))
        dsigdE_d = ((2.0*param.GF**2*m_n)/np.pi)*(p_d*(g_Ld**2+g_Rd**2*(Enu_2/Enu_1)**2)+p_d_b*(g_Rd**2+g_Ld**2*(Enu_2/Enu_1)**2))
        dsigdE = dsigdE_u + dsigdE_d
    elif param.neutype == "antineutrino":
        dsigdE_u = ((2.0*param.GF**2*m_n)/np.pi)*(p_u*(g_Ru**2+g_Lu**2*(Enu_2/Enu_1)**2)+p_u_b*(g_Lu**2+g_Ru**2*(Enu_2/Enu_1)**2))
        dsigdE_d = ((2.0*param.GF**2*m_n)/np.pi)*(p_d*(g_Rd**2+g_Ld**2*(Enu_2/Enu_1)**2)+p_d_b*(g_Ld**2+g_Rd**2*(Enu_2/Enu_1)**2))
        dsigdE = dsigdE_u + dsigdE_d        
    else :
        quit()
        
    return dsigdE[0]
    
def difsigDISCC(Enu_1,Enu_2,track,body,param):
    """ Calculates differential neutrino CC DIS cross section.
    
    d/dE_2_nu sig(\nu_1 N -> \nu_2 X)
    
    Ref.: Appendix B.  arXiv: hep-ph/0506298
    
    @type   Enu_1 	:	float
    @param  Enu_1	:	neutrino energy [GeV]
    @type   Enu_2 	:	float
    @param  Enu_2	:	neutrino energy [GeV]
    @type   track 	:	track object
    @param  track	:	trayectory
    @type   body 	:	body
    @param  body	:	body [Sun,Earth,Vacuum]
    @type   param       :	physics parameter set
    @param  param       :	physics parameter set

    @rtype	        :	float
    @return	        :	neutrino differential CC DIS cross section. [cm^2 eV^-1]
    """  
    return 0.0
    
#===============================================================================#
# Muon mean inelasticity                                                        #
#===============================================================================#

### Neutral and charge current mean inelasticity
## Ref : Arxiv hep-ph/9512364v1
# aprox. 1.0e0


def MuonMeanInelasticity_CC_Tbl(Enu,neu):
    """ Returns the muon charge current mean inelasticity from the tabulated values
    of the reference.
    
    Ref: arXiv hep-ph/9512364v1
    
    Note : value for 1GeV has been extrapolated.
        
    @type   Enu 	:	float
    @param  Enu 	:	neutrino energy [GeV]
    @type   neu 	:	int
    @param  neu 	:	0: neutrino, 1: antineutrino.

    @rtype	        :	float
    @return	        :	muon charge current mean inelasticity (y) [dimensionless]
    """  
    # extrapolating for 1GeV
    m = (0.483 - 0.477)/(1.0e2-1.0e1)
    ymuCC_1GeV = 0.487 + (1.0e0-1.0e1)*m

    if(neu == 0):
        E = [0.5,1.0e0,1.00e+001,1.00e+002,1.00e+003,1.00e+004,1.00e+005,1.00e+006,1.00e+007,1.00e+008,1.00e+009,1.00e+010,1.00e+011,1.00e+012]		
        yCC = [ymuCC_1GeV,ymuCC_1GeV,0.483,0.477,0.472,0.426,0.332,0.273,0.25,0.237,0.225,0.216,0.208,0.205]
        inter=interpolate.interp1d(E,yCC)
    elif(neu == 1):
        E = [0.0,1.0e0,1.00e+001,1.00e+002,1.00e+003,1.00e+004,1.00e+005,1.00e+006,1.00e+007,1.00e+008,1.00e+009,1.00e+010,1.00e+011,1.00e+012]	
        yCC = [0.0,0.0,0.333,0.340,0.354,0.345,0.301,0.266,0.249,0.237,0.225,0.216,0.208,0.205]				
        inter=interpolate.interp1d(E,yCC)
    else:
        print "NC:NEU:XSECTIONS:ERROR: MuonMeanInelasticity_CC_Tbl : Wrong neutrino type."
        quit()
    return inter(Enu)

def MuonMeanInelasticity_NC_Tbl(Enu,neu):
    """ Returns the muon neutral current mean inelasticity from the tabulated values
    of the reference.
    
    Ref: arXiv hep-ph/9512364v1
    
    Note : value for 1GeV has been extrapolated.
        
    @type   Enu 	:	float
    @param  Enu 	:	neutrino energy [GeV]
    @type   neu 	:	int
    @param  neu 	:	0: neutrino, 1: antineutrino.

    @rtype	        :	float
    @return	        :	muon neutral current mean inelasticity (y) [dimensionless]
    """  
    if(neu == 0):
        E = [1.00e+001,1.00e+002,1.00e+003,1.00e+004,1.00e+005,1.00e+006,1.00e+007,1.00e+008,1.00e+009,1.00e+010,1.00e+011,1.00e+012]		
        yNC = [0.474,0.470,0.467,0.428,0.341,0.279,0.254,0.239,0.227,0.217,0.210,0.207]
        inter=interpolate.interp1d(E,yNC)
    elif(neu == 1):
        E = [1.00e+001,1.00e+002,1.00e+003,1.00e+004,1.00e+005,1.00e+006,1.00e+007,1.00e+008,1.00e+009,1.00e+010,1.00e+011,1.00e+012]	
        yNC = [0.350,0.354,0.368,0.358,0.313,0.273,0.253,0.239,0.227,0.217,0.210,0.207]				
        inter=interpolate.interp1d(E,yNC)
    else:
        print "NC:NEU:XSECTIONS:ERROR: MuonMeanInelasticity_NC_Tbl : Wrong neutrino type."
        quit()
    return inter(Enu)
	

## EXPERIMENTAL

def MuonMeanInelasticity_CC_Avg(E,neutype,datapath = global_datapath):
  """ Mean charge-current muon inelasticity parameter averaged to 1.0 considering
  relationships between RES, QE and DIS cross sections.

  @type   E 	:	float
  @param  E	    :	neutrino energy [GeV]
  @type   neu	:	integer
  @param  neu	:	neu : 0 : neutrino, neu : 1 : antineutrino

  @rtype	:	float
  @return	:	charge current muon ineslasticity (y)	[dimensionless]
  """
  
  if E >= 100.0 or E < 1.0:
	return ymuCC(E,neutype)
  else : 
	existfile = False 
	for ff in os.listdir(datapath):
	  if ff == "y_avg_neutrino.dat" and neutype == 0:
	    existfile = True
	  elif ff == "y_avg_antineutrino.dat" and neutype == 1:
	    existfile = True

	if not existfile :       
	  QE_xsec_file  = "xsections_QE.dat"
	  RES_xsec_file = "xsections_RES.dat"
	  DIS_xsec_file = "xsections_DIS.dat"
	  
	  file_QE 	= open(datapath+QE_xsec_file,'r')
	  file_RES 	= open(datapath+RES_xsec_file,'r')
	  file_DIS 	= open(datapath+DIS_xsec_file,'r')
	  
	  h,QE  = gt.hreadfilev2(file_QE)
	  h,RES = gt.hreadfilev2(file_RES)
	  h,DIS = gt.hreadfilev2(file_DIS)
	  
	  QE_inter  = interpolate.interp1d(np.array(QE)[:,0],np.array(QE)[:,1])
	  RES_inter = interpolate.interp1d(np.array(RES)[:,0],np.array(RES)[:,1])
	  DIS_inter = interpolate.interp1d(np.array(DIS)[:,0],np.array(DIS)[:,1])
	  
	  Enu = np.arange(1.0,100.0,0.1) # [GeV]
	  ratio_DIS = map(lambda EE : DIS_inter(EE)/(DIS_inter(EE)+QE_inter(EE)+RES_inter(EE)), Enu)
	  ratio_QE_RES = map(lambda EE : (RES_inter(EE)+QE_inter(EE))/(DIS_inter(EE)+QE_inter(EE)+RES_inter(EE)), Enu)
	  
	  avg_y = []
	  avg_y_list =[]
	  for i,EE in enumerate(Enu):
	    avg_y.append(ymuCC(EE,neutype)*ratio_DIS[i]+0.0*ratio_QE_RES[i])
	    avg_y_list.append([EE,ymuCC(EE,neutype)*ratio_DIS[i]+0.0*ratio_QE_RES[i]])
	    
	  if neutype == 0:
	    file = open(datapath+"y_avg_neutrino.dat", 'w')
	  elif neutype == 1:
	    file = open(datapath+"y_avg_antineutrino.dat", 'w')
	  gt.quickprint(file,avg_y_list)
	else :
	  if neutype == 0:
	    file = open(datapath+"y_avg_neutrino.dat", 'r')
	  elif neutype == 1:
	    file = open(datapath+"y_avg_antineutrino.dat", 'r')
	  h,dat = gt.hreadfilev2(file)
	  
	  dat = np.array(dat)
	  Enu = dat[:,0]
	  avg_y = dat[:,1]
	if E > Enu[-1]:
	  return ymuCC(E,neutype)
	else :
	  inter_y_avg = interpolate.interp1d(Enu,avg_y)
	  return inter_y_avg(E)
	  
#===============================================================================#
# MINOS CROSS SECTIONS                                                          #
#===============================================================================#

def nuMINOSxsection_CC(E,neu):
  """ neutrino-nucleon inclusive CC cross section
  
  # REFERENCE ARXIV 0910.2201v2 pag. 27
  @type   E 	:	float
  @param  E	:	neutrino energy [GeV]
  @type   neu	:	integer
  @param  neu	:	neu : 0 : neutrino, neu : 1 : antineutrino

  @rtype	:	float
  @return	:	crosssection nu-N. New : [eV^-2] Old : [cm^2].
  """
  
  Ev 		= [1.0,2.0,3.48,4.45,5.89,7.97,10.45,13.43,16.42,19.87,23.88,27.87,32.82,38.87,45.77]
  Eav 		= [6.07,7.99,10.43,13.42,16.41,19.82,23.82,27.84,32.72,38.74,45.61]
  signeu 	= [0.925,0.8,0.748,0.711,0.708,0.722,0.699,0.691,0.708,0.689,0.683,0.686,0.675,0.675,0.676]
  siganeu 	= [0.305,0.300,0.303,0.314,0.304,0.316,0.320,0.332,0.325,0.352,0.324]
  if neu == 0 :
    inter=interpolate.interp1d(Ev,signeu)
    if E<1.0:
      return 0.925*E*1.0e-38*pc.cm**2
    elif E>45.77:
      return 0.676*E*1.0e-38*pc.cm**2
    else:
      return inter(E)*E*1.0e-38*pc.cm**2
  elif neu ==1 :
    inter=interpolate.interp1d(Eav,siganeu)
    if E<6.07:
      return 0.305*E*1.0e-38*pc.cm**2
    elif E>45.61:
      return 0.324*E*1.0e-38*pc.cm**2
    else:
      return inter(E)*E*1.0e-38*pc.cm**2
  else :
    print "Wrong neutrino type."
    quit()

def nuMINOSxsection_CC_binned(E,neu):
  """ neutrino-nucleon inclusive CC cross section
  
  # REFERENCE ARXIV 0910.2201v2 pag. 27
  @type   E 	:	float
  @param  E	:	neutrino energy [GeV]
  @type   neu	:	integer
  @param  neu	:	neu : 0 : neutrino, neu : 1 : antineutrino

  @rtype	:	float
  @return	:	crosssection nu-N [cm^2]
  # sig/E [10^{-38}cm^2/GeV]
  """
  if neu == 0 :
    if E < 3.0 :
      sig = 0.0
    elif E>=3.0 and E<=4.0 :
      sig = 0.748
    elif E>4.0 and E<=5.0 :
      sig = 0.711
    elif E>5.0 and E<=7.0 :
      sig = 0.708
    elif E>7.0 and E<=9.0 :
      sig = 0.722
    elif E>9.0 and E<=12.0 :
      sig = 0.699
    elif E>12.0 and E<=15.0 :
      sig = 0.691
    elif E>15.0 and E<=18.0 :
      sig = 0.708
    elif E>18.0 and E<=22.0 :
      sig = 0.689
    elif E>22.0 and E<=26.0 :
      sig = 0.683
    elif E>26.0 and E<=30.0 :
      sig = 0.686
    elif E>30.0 and E<=36.0 :
      sig = 0.675
    elif E>36.0 and E<=42.0 :
      sig = 0.675
    elif E>42.0 and E<=50.0 :
      sig = 0.675
    elif E>50.0:
      sig = 0.675
  elif neu == 1:
    if E <= 3.0 :
      sig = 0.0
    elif E>3.0 and E<=4.0 :
      sig = 0.305#cero
    elif E>4.0 and E<=5.0 :
      sig = 0.305#cero
    elif E>5.0 and E<=7.0 :
      sig = 0.305
    elif E>7.0 and E<=9.0 :
      sig = 0.300
    elif E>9.0 and E<=12.0 :
      sig = 0.303
    elif E>12.0 and E<=15.0 :
      sig = 0.314
    elif E>15.0 and E<=18.0 :
      sig = 0.304
    elif E>18.0 and E<=22.0 :
      sig = 0.316
    elif E>22.0 and E<=26.0 :
      sig = 0.320
    elif E>26.0 and E<=30.0 :
      sig = 0.332
    elif E>30.0 and E<=36.0 :
      sig = 0.325
    elif E>36.0 and E<=42.0 :
      sig = 0.352
    elif E>42.0 and E<=50.0 :
      sig = 0.324
    elif E>50.0:
      sig = 0.324
  else:
    print "Wrong neutrino type."
    quit()
  return 1.0e-38*sig*E
  
#===============================================================================
# MiniBooNE cross sections
#===============================================================================

def LoadMiniBooNECS(mbparam):
    """ Reads MiniBooNE cross sections and stores it in memory.
    
    NOT CURRENTLY BEENING USED. MISSING IMPLEMENTATION.

    @type   mbparam   :   miniboone_config
    @param  mbparam   :   miniboone run configuration
    
    @rtype            :   array
    @return           :   raw cross section data array
    """
    if mbparam.CrossSection.LoadCrossSection == []:
      file  = open(datapath + "CC_charged_pion_XSecTables.dat",'r')
      h,dat = gt.hreadfilev3(file)
      file.close()
      mbparam.CrossSection.LoadCrossSection = dat
      return dat
    else :
      return mbparam.CrossSection.LoadCrossSection
      
def MiniBooNECS(E,mbparam):
    """ Returns MiniBooNE CS for a given detection channel.
    
      - Charged current neutral pion CS
      Ref : A.A. Aguilar-Arevalo et al., "Measurement of Muon Neutrino Induced Charged Current Neutral Pion Cross Sections on Mineral Oil at Enu=0.5-2.0 GeV"
      arXiv:1010.3264 [hep-ex]
      
      - Neutral current neutral pion CS
      Ref : A.A. Aguilar-Arevalo et al., "Measurement of νμ and ν̅μ induced neutral current single π0 production cross sections on mineral oil at Eν~O(1 GeV)", arXiv:0911.2063 [hep-ex]A.A. Aguilar-Arevalo et al., "Measurement of νμ and ν̅μ induced neutral current single π0 production cross sections on mineral oil at Eν~O(1 GeV)"
      arXiv:0911.2063 [hep-ex]
      
      - Charged current charged pion CS
      Ref : A.A. Aguilar-Arevalo et al., "Measurement of Neutrino-Induced Charged Current-Charged Pion Production Cross Sections on Mineral Oil at Enu ~ 1 GeV"
      arXiv:1011.3572 [hep-ex]    
      
      - Charged current quasielastic CS
      Ref : A.A. Aguilar-Arevalo et al., "First Measurement of the Muon Neutrino Charged Current Quasielastic Double Differential Cross section"
      arXiv:1002:2680 [hep-ex]

    @type   E         :   float
    @param  E         :   (anti)neutrino energy         [eV]
    @type   mbparam   :   miniboone_config
    @param  mbparam   :   miniboone run configuration
    
    @rtype            :   array
    @return           :   cross section data array [cm^2/nucleon]
    """
    
    if mbparam.neutype == "neutrino":
      ## Energies in GeV. Cross Section in 1.0e-39 cm^2 / nucleon
      # CH_2
      neutron_number_CH2 = 6.0
      proton_number_CH2  = 6.0+2.0
      nucleon_number_CH2 = 6.0+6.0+2.0
      
      # some of the cross section where not measured per nucleon, but per CH_2.
      # to convert to per nucleon cross section we rescale them appropiartly
      # beware that this rescaling does not take account of nuclear effects
      
      ## CC_neutral_pion
      Enu_0           = [0.50,0.60,0.70,0.80,0.90,1.00,1.10,1.20,1.30,1.40,1.50,1.60,1.70,1.90]
      CC_neutral_pion = [1.76,3.83,5.68,7.31,9.20,11.06,12.42,13.89,15.23,16.38,18.20,19.37,20.80,21.92]
      CC_neutral_pion_nucleon = map(lambda x : x*(neutron_number_CH2/nucleon_number_CH2)/nucleon_number_CH2,CC_neutral_pion)
      ## NC_neutral_pion
      Enu_1           = [0.0218023255814, 0.0726744186047, 0.15261627907, 0.232558139535, 0.297965116279, 0.348837209302, 0.392441860465, 0.428779069767, 0.479651162791, 		0.508720930233, 0.545058139535, 0.588662790698, 0.654069767442, 0.704941860465, 0.755813953488, 0.828488372093, 0.893895348837, 0.952034883721, 0.988372093023, 1.06104651163,  1.13372093023, 1.18459302326, 1.25, 1.31540697674, 1.36627906977, 1.41715116279, 1.46802325581, 1.51889534884, 1.5625, 		1.6351744186, 1.69331395349, 1.74418604651, 1.80959302326, 1.85319767442, 1.88226744186, 1.93313953488, 1.98401162791, 2.04215116279, 2.09302325581, 		2.13662790698, 2.20930232558, 2.27470930233, 2.35465116279, 2.40552325581, 2.4636627907]
      NC_neutral_pion = [0.0353975804557, 0.0304498473103, 0.0356618510688, 0.0358233497768, 0.0359554850834, 0.041108762039, 0.0512978623444, 0.07157329105, 0.107029598309, 0.127290345314, 0.167767794221, 0.218360934931, 0.304351656096, 0.385262508809, 0.451021846371, 0.542077754287, 0.623017970402, 	0.688791989664, 0.734319943622, 0.800123326286, 0.876027719051, 0.921585036411, 0.977272727273, 1.02790991308, 1.06336622034, 1.0988225276, 1.13427883486, 	1.16468463707, 1.19002525253,  1.22552560489, 1.25089558374, 1.2762508809, 1.30668604651, 1.33202666197, 1.34723690392, 1.36249119098, 1.38784648814,	1.40816596195, 1.42847075405, 1.43865985436, 1.46910970167, 1.48439335213, 1.50980737609, 1.52001115809, 1.53528012685]
      NC_neutral_pion_nucleon = map(lambda x : x/nucleon_number_CH2,NC_neutral_pion)
      ## CC_charged_pion
      Enu_2           = [0.50,0.60,0.65,0.70,0.75,0.80,0.85,0.90,0.95,1.00,1.05,1.10,1.15,1.20,1.25,1.30,1.35,1.40,1.45,1.50,1.55,1.60,1.65,1.70,1.75,1.80,1.90]
      CC_charged_pion = [6.1,11.5,15.8,19.6,24.1,28.3,32.5,37.3,41.6,46.6,49.7,52.9,56.3,59.1,62.3,66.7,70.3,72.3,77.6,80.8,83.7,86.4,88.5,93.0,92.6,97.1,99.2]
      CC_charged_pion_nucleon = map(lambda x : x/nucleon_number_CH2,CC_charged_pion)
      ## CC_QE
      Enu_3           = [0.40,0.45,0.50,0.55,0.60,0.65,0.70,0.75,0.80,0.90,1.0,1.1,1.3,1.5,2.0]
      CC_QE           = [7.985,8.261,8.809,9.530,10.13,10.71,11.11,11.55,12.02,12.30,12.58,12.58,12.78,12.36]
      CC_QE_nucleon           = map(lambda x : x*(neutron_number_CH2/nucleon_number_CH2),CC_QE)

    elif mbparam.neutype == "antineutrino":
      ## Energies in GeV. Cross Section in 1.0e-39 cm^2 / nucleon
      ## CC_neutral_pion
#      Enu_0           = []
#      CC_neutral_pion = []
      ## NC_neutral_pion
      Enu_1           = [0.315942028986,0.344918155809,0.36808780076,0.36808780076,0.36808780076,0.399934337039,0.423085221143,0.443337554524,0.457802166878,0.469349467661,0.483814080015,0.506918062005,0.53580976502,0.56468270719,0.590666479058,0.622409830683,0.63973547207,0.65704235261,0.668570892547,0.685887153511,0.708962994231,0.734937385676,0.763791567,0.781107827963,0.795516157779,0.815684067351,0.84166783922,0.861873270484,0.879161390179,0.902255991745,0.919553491862,0.936869752826,0.95415787252,0.971474133483,0.983002673421,0.997448524928]
      NC_neutral_pion = [0,0.00970873786408,0.0291262135922,0.0291262135922,0.0291262135922,0.0679611650485,0.106796116505,0.145631067961,0.174757281553,0.223300970874,0.252427184466,0.339805825243,0.436893203883,0.553398058252,0.660194174757,0.805825242718,0.873786407767,0.961165048544,1.02912621359,1.1067961165,1.22330097087,1.33980582524,1.47572815534,1.55339805825,1.64077669903,1.76699029126,1.87378640777,1.96116504854,2.06796116505,2.16504854369,2.26213592233,2.33980582524,2.44660194175,2.52427184466,2.59223300971,2.64077669903]
      ## CC_charged_pion
#      Enu_2           = []
#      CC_charged_pion = []
      ## CC_QE
      Enu_3           = [0.131947552168,0.135031403787,0.138187330565,0.141417017029,0.148104604888,0.158733611101,0.170125427985,0.186596292318,0.199987727054,0.214340223351,0.235091787546,0.296191032977,0.373169683773,0.470154722402,0.60618989935,0.729249691684,0.918778242886,1.08005237452,1.32967596747,1.52737785039,1.71440613585,1.96931133791,2.26211693047,2.53911443368,2.85003043852,3.27378509595]
      CC_QE           = [0.118230563003,0.126675603217,0.137533512064,0.148391420912,0.160455764075,0.174932975871,0.185790884718,0.196648793566,0.207506702413,0.212332439678,0.219571045576,0.224396782842,0.223190348525,0.220777479893,0.215951742627,0.207506702413,0.200268096515,0.190616621984,0.17855227882,0.167694369973,0.158042895442,0.147184986595,0.136327077748,0.125469168901,0.117024128686]
    else :
      print "Wrong neutrino type."
      quit()
    
    # Calculating all the CS for this energy
    if mbparam.CrossSection.CalculationMode == "binned":
      cs_0 = gt.BinnedValue(E/mbparam.GeV,Enu_0,CC_neutral_pion_nucleon)
      cs_1 = gt.BinnedValue(E/mbparam.GeV,Enu_1,NC_neutral_pion_nucleon)
      cs_2 = gt.BinnedValue(E/mbparam.GeV,Enu_2,CC_charged_pion_nucleon)
      cs_3 = gt.BinnedValue(E/mbparam.GeV,Enu_3,CC_QE_nucleon)
      return [cs_0*1.0e-39,cs_1*1.0e-39,cs_2*1.0e-39,cs_3*1.0e-39]
    elif mbparam.CrossSection.CalculationMode == "interpolated":
      inter_0 = interpolate.interp1d(Enu_0,CC_neutral_pion_nucleon)
      inter_1 = interpolate.interp1d(Enu_1,NC_neutral_pion_nucleon)
      inter_2 = interpolate.interp1d(Enu_2,CC_charged_pion_nucleon)
      inter_3 = interpolate.interp1d(Enu_3,CC_QE_nucleon)
      return [inter_0(E/mbparam.GeV)*1.0e-39,inter_1(E/mbparam.GeV)*1.0e-39,inter_2(E/mbparam.GeV)*1.0e-39,inter_3(E/mbparam.GeV)*1.0e-39]
    else:
      print "Wrong cross section calculation mode."
      quit()
  
	  

if __name__ == '__main__':
    pass
