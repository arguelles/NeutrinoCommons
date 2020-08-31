from numpy import log as log
# python modules
import os as os
# my modules
import DM as DM
import numpy as np
import generaltools as gt
import physicsconstants as PC

## external C functions ##

#cdef extern from "math.h":
#    double log(double)



## main functions ##

cpdef DMFNeuFluxMCDetOpt(ch,DMm,DMsig,p,oo,dp,cc):
    #,param,onlyosc,datapath,crosscheck):
    param = p
    onlyosc = oo
    datapath = dp
    crosscheck = cc
    
    # variable definition
    cdef double DM_annihilation_rate_Sun, normalization
    cdef double DMm_GeV
    
    cdef int nu_bin_num
    cdef double Emin,Emax
    cdef double point_num
    cdef double x
    
    cdef list[double] E_nu_list,E_bin_width,E_nu_hpl,E_nu_bin,E_anu_bin,E_bin_ratio
    
    #cdef list[string] files
    #cdef string filename,MCdatapath
    
    cdef int ineu,i,neuneu
    cdef int family,ii,jj,iineu,ini
    
    cdef double ini_evt_energy,E_nu_in,E_nu_out,E_nu_list_0,log_E_bin_ratio

    ## BEGIN CREATING BINS ##
    DM_annihilation_rate_Sun = DM.DMSunAnnihilationRate(DMm,DMsig,param)           # [eV]
    normalization = np.sum((DM.DM_annihilation_rate_Sun/(4.0*np.pi*param.AU**2)))  # [eV^3]
    
    DMm_GeV = DMm/param.GeV
    
    #print normalization
    #quit()
    # assuming neutrino binnum = 30
    nu_bin_num  = 200
    point_num   = 1000.0
    Emin        = 1.0
    Emax        = 10000.0
    
    E_nu_list   = gt.LogSpaceEnergies(Emin,Emax,binnum = nu_bin_num)
    E_bin_width = [E_nu_list[i+1]-E_nu_list[i] for i in range(len(E_nu_list)-1)]
    E_nu_hpl    = gt.MidPoint(gt.LogSpaceEnergies(Emin,Emax,binnum = nu_bin_num))    
    E_nu_bin    = [0.0]*nu_bin_num # neutrino bins
    E_anu_bin   = [0.0]*nu_bin_num # antineutrino bins
    E_bin_ratio = E_nu_list[1]/E_nu_list[0]
    
    log_E_bin_ratio = log(E_bin_ratio)
    E_nu_list_0 = E_nu_list[0]    
    
    ## END CREATING BINS ##
    
    DM_pdf = []
    for neutype in range(6):
        DM_pdf.append(DM.DM_distribution(ch,DMm/param.GeV,neutype))
            
    # calculate DM distribution arrays for futher use
    DM_pdf_table = []
    for neuneu in range(6):
        DM_pdf_table.append(map(lambda EE : float(DM_pdf[neuneu].PDF(EE)),E_nu_hpl))    
    
    for ineu in [2]:#range(3):
        ## BEGIN READING DATA FROM MC ## 
        # changed for no osc
        if onlyosc :
            MCdatapath = datapath+"legion_ineu_"+str(ineu)+"_"+param.name+"_noosc/"
        else :
            MCdatapath = datapath+"legion_ineu_"+str(ineu)+"_"+param.name+"/"
        rparam  = PC.PhysicsConstants()
        
        files = []
        for filename in os.listdir(MCdatapath):
            files.append(filename)
                    
        # load all events
        evt = []
        for filename in files :
            print MCdatapath
            print filename
            file = open(MCdatapath+filename,'r')
            data = []
            gt.hreadfilev4(file,data,rparam)
            if crosscheck :
                if gt.Compareparams(param,rparam):
                    print "Using : "+filename
                    for e in data :
                        for ee in e:
                            evt.append(ee)
            else :
                print "Using : "+filename
                for e in data :
                    for ee in e:
                        evt.append(ee)                        
        
        ## END READING DATA FROM MC ##
        
        # GET DARK MATTER DISTRIBUTION
        #flavor of the neutrino (0 : nu_e, 1 : anu_e, 2 : nu_mu 0, 3 : anu_mu, 4 : nu_tau, 5 : anu_tau)

           
        for i,e in enumerate(evt):
            if i > 1000 :
                break
            #print i,len(evt),str((i/len(evt))*100) + " %"
            if len(e) > 5:
                
                if param.name == "STD":
                    
                    neutrino = True
                    
                    family = e[0]
                    #try:
                    #    next_family = evt[i+1]
                    #    # cross check this
                    #    if family == next_family and e[1] != 2 :
                    #        neutrino = False
                    #except:
                    #    pass
                    #
                    #try:
                    #    ii = 1
                    #    same_evt = True
                    #    while same_evt :
                    #        next_family = evt[i+ii]
                    #        if len(next_family) == 3 and evt[i+ii][0] != evt[i+ii+1]:
                    #            break
                    #        else :
                    #            ii = ii + 1
                    #            
                    #    ii = ii-1
                    #except:
                    #    ii = 0
                    
                    ini_evt_energy = evt[i+ii][2]
                    
                    iineu = e[1]
                    E_nu_in  = e[2]
                    E_nu_out = e[3]
                    
                    jj = int(log(E_nu_in/E_nu_list_0)/log_E_bin_ratio)
                    ii = int(log(E_nu_out/E_nu_list_0)/log_E_bin_ratio)
                    
                    ini = int(log(ini_evt_energy/E_nu_list_0)/log_E_bin_ratio)
                    
                    if neutrino:
                        DMpdf = DM_pdf[iineu*2]
                        if E_nu_in + E_bin_width[jj] > DMm_GeV :
                            E_nu_bin[ii]  = E_nu_bin[ii]  + e[5]*(DM_pdf_array[ini]/DMm_GeV)*(DMm_GeV-E_nu_list[ini])
                        else:
                            E_nu_bin[ii]  = E_nu_bin[ii]  + e[5]*(DM_pdf_array[ini]/DMm_GeV)*E_bin_width[ini]
                    else :
                        DMpdf = DM_pdf[iineu*2+1]
                        E_anu_bin[ii] = E_anu_bin[ii] + e[5]*(DM_pdf_array[ini]/DMm_GeV)*E_bin_width[ini]
                        
                else :
                    family = e[0]
                    
                    #try:
                    #    ii = 1
                    #    same_evt = True
                    #    while same_evt :
                    #        next_family = evt[i+ii]
                    #        if len(next_family) == 4 and evt[i+ii][0] != evt[i+ii+1]:
                    #            break
                    #        else :
                    #            ii = ii + 1
                    #            
                    #    ii = ii-1
                    #except:
                    #    ii = 0
                    
                    ini_evt_energy = evt[i+ii][3]
                    
                    iineu = e[1]
                    ineutype = e[2]
                    E_nu_in  = e[3]
                    E_nu_out = e[4]
                    
                    log_E_bin_ratio = np.log(E_bin_ratio)
                    E_nu_list_0 = E_nu_list[0]
                    
                    jj = int(np.log(E_nu_in/E_nu_list_0)/log_E_bin_ratio)
                    ii = int(np.log(E_nu_out/E_nu_list_0)/log_E_bin_ratio)
                    
                    ini = int(np.log(ini_evt_energy/E_nu_list_0)/log_E_bin_ratio)
                    
                    if ineutype == 0:
                        DM_pdf_array = DM_pdf_table[iineu*2]
                        if E_nu_in + E_bin_width[jj] > DMm_GeV :
                            E_nu_bin[ii]  = E_nu_bin[ii]  + e[6]*(DM_pdf_array[ini]/DMm_GeV)*(DMm_GeV-E_nu_list[ini])
                        else:
                            E_nu_bin[ii]  = E_nu_bin[ii]  + e[6]*(DM_pdf_array[ini]/DMm_GeV)*E_bin_width[ini]
                    elif ineutype == 1 :
                        DM_pdf_array = DM_pdf_table[iineu*2+1]
                        if E_nu_in + E_bin_width[jj] > DMm_GeV :
                            E_anu_bin[ii] = E_anu_bin[ii] + e[6]*(DM_pdf_array[ini]/DMm_GeV)*(DMm_GeV-E_nu_list[ini])
                        else :
                            E_anu_bin[ii] = E_anu_bin[ii] + e[6]*(DM_pdf_array[ini]/DMm_GeV)*E_bin_width[ini]

    E_nu_bin = [normalization*x/(point_num) for x in E_nu_bin]
    E_anu_bin = [normalization*x/(point_num) for x in E_anu_bin]    
    
    return [E_nu_hpl,E_nu_bin,E_anu_bin]  