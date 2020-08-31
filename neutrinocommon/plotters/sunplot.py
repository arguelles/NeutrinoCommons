import numpy as np
import zipfile as zp
import os
import re
import fnmatch
import scipy.interpolate as interpolate
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.offsetbox as osb
import matplotlib.patches as ptc
import pylab as pl
# my modules
import neutrinocommon.astro.body as bd
import neutrinocommon.astro.DM
import neutrinocommon.neu.neuosc as no
import neutrinocommon.neu.xsections as xs
import neutrinocommon.physconst.physicsconstants as PC
import neutrinocommon.tools.generaltools as gt
import neutrinocommon.exp.icecube as ice

#===============================================================================
# SUN MATTER PROFILE
#===============================================================================

def PlotSolarParameters():
    plt.figure()
    R = np.arange(0.01,1.0,0.01)
    Sun = bd.Sun()
    
    ye = map(lambda r : float(Sun.rye(r)),R)
    xh = map(lambda r : float(Sun.rxh(r)),R)
    negrnc = map(lambda r : -float(Sun.rRNC(r)),R)
    
    pye = plt.plot(R,ye,color = 'black',linestyle = 'dashed')
    pxh = plt.plot(R,xh,color = 'black',linestyle = 'solid')
    pnegrnc = plt.plot(R,negrnc,color = 'black',linestyle = 'dotted')
    plt.legend([pye,pxh,pnegrnc],[r"$y_e$",r"$X_h$","$-R_{NC}$"],loc = 7)
    plt.axis([0.0, 1.0, 0.0, 1.0])
    plt.savefig("PlotSolarParameters.png")
    
#===============================================================================
# STERILE NEUTRINO OSCILLATION PROBABILITY
#===============================================================================

def PlotNeuOscProb(ineu,fneu,Ri,Enumin,Enumax,body,param,rkcomp):
    """Plots P(neu_ineu -> neu_fneu) as a function of the Energy
    # ineu,fneu     : 0 (electron), 1 (muon), 2 (tau)
    # Ri            : neutrino production point in the body
    # Enumin        : minimum neutrino energy
    # Enumax        : maximum neutrino energy
    # body          : body where the neutrino propagates
    # param         : physics parameter set
    # rkcomp        : toggle RK validation
    # If rkcomp = True, will use RK to compute the probability and compare
    """
    plt.figure()
    ordmag = np.log10(Enumax)-np.log10(Enumin)
    npoints = 1000.0*ordmag
    # Figuring out best energy scale
    try:
        if(Enumax/param[0].MeV <= 500.0) :
            scale =  param[0].MeV
            scalen = "MeV"
        elif(Enumax/param[0].GeV <= 500.0) :
            scale =  param[0].GeV
            scalen = "GeV"
        else :
            scale =  param[0].TeV
            scalen = "TeV"
    except (TypeError,AttributeError):
        if(Enumax/param.MeV <= 500.0) :
            scale =  param.MeV
            scalen = "MeV"
        elif(Enumax/param.GeV <= 500.0) :
            scale =  param.GeV
            scalen = "GeV"
        else :
            scale =  param.TeV
            scalen = "TeV"
    # Figuring out ylabel
    if (ineu == 0):
        ineulabel = "e"
        sineulabel = "e"
    elif (ineu == 1):
        ineulabel = "\mu"
        sineulabel = "mu"
    elif (ineu == 2):
        ineulabel = "\\tau"
        sineulabel = "tau"
    elif (ineu > 2) :
        ineulabel = "s"
        sineulabel = "s"
    if (fneu == 0):
        fneulabel = "e"
        sfneulabel = "e"
    elif (fneu == 1):
        fneulabel = "\mu"
        sfneulabel = "mu"
    elif (fneu == 2):
        fneulabel = "\\tau"
        sfneulabel = "tau"
    elif (fneu > 2) :
        fneulabel = "s"
        sfneulabel = "s"
            
    # RK points
    ERKstep = (np.log10(Enumax)-np.log10(Enumin))/(10.0)    
    ERK = np.arange(np.log10(Enumin),np.log10(Enumax),ERKstep)
    ERK = map(lambda E : (10**E)/scale,ERK)
    
    Estep = (Enumax-Enumin)/npoints        
    Enu = np.arange(Enumin,Enumax,Estep)/scale
    # Generating plot
    try:
        pid = [0]*len(param)
        for i,p in enumerate(param):
            p.Refresh()
            PMNS = no.mixmatrix(p)
            fM2 = no.flavorM2(p)
            Rf = body.Radius*p.km
            track = body.track(Ri,Ri,Rf)
            Pe = map(lambda E : no.AdiabaticProbability(ineu,fneu,E*scale,Ri,body,PMNS,fM2,p),Enu)
            plt.xlabel(r"$\mathrm{E}_\nu\mathrm{["+scalen+"]}$")
            ylabel = "$ P(\\nu_{"+ineulabel+"} \\rightarrow \\nu_{"+fneulabel+"}) $"
            plt.ylabel(r""+ylabel)
            pid[i]= plt.plot(Enu,Pe,color = 'black',label = r""+p.name,linestyle = p.linestyle)
            if rkcomp :
                PeRK = map(lambda E : no.AvgNeuProb_RK_STD(ineu,fneu,E*scale,p),ERK)
                #PeRK = map(lambda E : no.probneu(ineu,fneu,E*scale,track,body,fM2,p),ERK)
                for i,e in enumerate(ERK):
                    plt.plot([e],[PeRK[i]],p.markerstyle,color='red',markersize=12,markeredgecolor='red',markeredgewidth=1,label='_nolegend_')
    except (TypeError,AttributeError):
        param.Refresh()
        PMNS = no.mixmatrix(param)
        fM2 = no.flavorM2(param)
        Rf = body.Radius*param.km
        track = body.track(Ri,Ri,Rf)
        Pe = map(lambda E : no.AdiabaticProbability(ineu,fneu,E*scale,Ri,body,PMNS,fM2,param),Enu)
        plt.xlabel(r"$\mathrm{E}_\nu\mathrm{["+scalen+"]}$")
        ylabel = "$ P(\\nu_{"+ineulabel+"} \\rightarrow \\nu_{"+fneulabel+"}) $"
        plt.ylabel(r""+ylabel)
        plt.plot(Enu,Pe,color = 'black',linestyle = param.linestyle)
        if rkcomp :
            PeRK = map(lambda E : no.probneu(ineu,fneu,E*scale,track,body,fM2,param),ERK)
            for i,e in enumerate(ERK):
                plt.plot([e],[PeRK[i]],param.markerstyle,color='red',markersize=12,markeredgecolor='red',markeredgewidth=1,label='_nolegend_')
    plt.legend(loc = 0)
    plt.semilogx()
    plt.axis([Enumin/scale,Enumax/scale,0.0,1.0])
    path = "../plots/"
    filename = path+"PlotNeuOscProb_Emin_"+str(Enumin/scale)+"_"+scalen+"_Emax_"+str(Enumax/scale)+"_"+scalen+"_nui_"+sineulabel+"_nuf_"+sfneulabel+".png"
    plt.savefig(filename)
    
    
def PlotSunAbsorptionProbability(Enumin,Enumax,param,datapath = "../data/NeutrinoAbsorptionProb/"):
    """Plots absorption probability due to CC interactions using analitic calculation.
    # iineu         : 0 (electron), 1 (muon), 2 (tau)
    # Enumin        : minimum neutrino energy       [eV]
    # Enumax        : maximum neutrino energy       [eV]
    # param         : physics parameter
    """
    
    E = gt.MidPoint(gt.LogSpaceEnergies(Enumin,Enumax,binnum = 200))
    
    filename_neutrino       = "DatNeuSunAbsorptionProbability_1.0_GeV_1000000.0_GeV_neutrino.dat"
    filename_antineutrino   = "DatNeuSunAbsorptionProbability_1.0_GeV_1000000.0_GeV_antineutrino.dat"
    
    file_neu  = open(datapath+filename_neutrino,'r')
    file_aneu = open(datapath+filename_antineutrino,'r')
    
    p_neu = []
    p_aneu = []   
    
    gt.hreadfilev4(file_neu,p_neu,param)
    gt.hreadfilev4(file_aneu,p_aneu,param)
    
    file_neu.close()
    file_aneu.close()
    
    E = map(lambda x : x[0], p_neu[0])
    p_neu = map(lambda x : x[1], p_neu[0])
    p_aneu = map(lambda x : x[1], p_aneu[0])
    
    #print p_neu
    
    fig = plt.figure()
    ax = plt.subplot(111)
    
    mpl.rcParams['axes.labelsize'] = "large"
    mpl.rcParams['xtick.labelsize'] = "large"
    mpl.rcParams['legend.fontsize'] = "small"
    
    mpl.rcParams['font.size'] = 15
    
    E_TeV = map(lambda x : x/param.TeV,E)
    
    plt.plot(E_TeV,p_neu,label = r"$\nu$", lw = 4)
    plt.plot(E_TeV,p_aneu,label = r"$\bar{\nu}$", lw = 4)
    
    plt.xlabel(r"$E_\nu \, [\mathrm{TeV}]$")
    plt.ylabel(r"$P_{absorption}$")
    
    plt.loglog()
    
    plt.legend(loc = "upper right")
    
    path = "../plots/"
    
    filename = "SunAbsorptionProbability.eps"
    
    plt.savefig(path+filename)
    
def PlotSunAbsorptionProbabilityMC(iineu,Enumin,Enumax,param,sparam = None,datapath = "../data/SunOscProbabilities/"):
    """Plots absorption probability due to CC interactions using MC data.
    # Enumin        : minimum neutrino energy       [eV]
    # Enumax        : maximum neutrino energy       [eV]
    # param         : physics parameter
    """
    
    E = gt.MidPoint(gt.LogSpaceEnergies(Enumin,Enumax,binnum = 200))
    neu_label  = ["$\\nu_e$","$\\nu_\\mu$","$\\nu_\\tau$"]
    aneu_label = ["$\\bar{\\nu}_e$","$\\bar{\\nu}_\\mu$","$\\bar{\\nu}_\\tau$"]
    
    ## setup canvas ##
    
    fig = plt.figure()
    ax = plt.subplot(111)
            
    mpl.rcParams['axes.labelsize'] = "large"
    mpl.rcParams['xtick.labelsize'] = "large"
    mpl.rcParams['legend.fontsize'] = "small"
            
    mpl.rcParams['font.size'] = 15
    
    ## end setup canvas ##
    
    ## FOR PARAM ##
     
    for ineu in [iineu]:
        
        filename_neutrino       = "DataNeuOscProb_RK_neutrino_Emin_"+str(Enumin/param.GeV)+"_GeV_Emax_"+str(Enumax/param.GeV)+"_GeV_ineu_"+str(ineu)+"_param_"+param.name+".dat"
        filename_antineutrino   = "DataNeuOscProb_RK_antineutrino_Emin_"+str(Enumin/param.GeV)+"_GeV_Emax_"+str(Enumax/param.GeV)+"_GeV_ineu_"+str(ineu)+"_param_"+param.name+".dat"
        
        file_neu  = open(datapath+filename_neutrino,'r')
        file_aneu = open(datapath+filename_antineutrino,'r')
        
        p_neu = []
        p_aneu = []   
        
        gt.hreadfilev4(file_neu,p_neu,param)
        gt.hreadfilev4(file_aneu,p_aneu,param)
        
        file_neu.close()
        file_aneu.close()
        
        E = map(lambda x : x[0], p_neu[0])
        p_neu  = map(lambda x : x[param.numneu + 1], p_neu[0])
        p_aneu = map(lambda x : x[param.numneu + 1], p_aneu[0])
        
        E_TeV = map(lambda x : x/param.TeV,E)
        
        plt.plot(E_TeV,p_neu,label = neu_label[ineu], lw = 4, linestyle = "solid")
        plt.plot(E_TeV,p_aneu,label = aneu_label[ineu], lw = 4, linestyle = "solid")
        
    ## FOR SPARAM ##
    if sparam != None :
        for ineu in [iineu]:
            
            filename_neutrino       = "DataNeuOscProb_RK_neutrino_Emin_"+str(Enumin/sparam.GeV)+"_GeV_Emax_"+str(Enumax/sparam.GeV)+"_GeV_ineu_"+str(ineu)+"_param_"+sparam.name+".dat"
            filename_antineutrino   = "DataNeuOscProb_RK_antineutrino_Emin_"+str(Enumin/sparam.GeV)+"_GeV_Emax_"+str(Enumax/sparam.GeV)+"_GeV_ineu_"+str(ineu)+"_param_"+sparam.name+".dat"
            
            file_neu  = open(datapath+filename_neutrino,'r')
            file_aneu = open(datapath+filename_antineutrino,'r')
            
            p_neu = []
            p_aneu = []   
            
            gt.hreadfilev4(file_neu,p_neu,sparam)
            gt.hreadfilev4(file_aneu,p_aneu,sparam)
            
            file_neu.close()
            file_aneu.close()
            
            E = map(lambda x : x[0], p_neu[0])
            p_neu  = map(lambda x : x[sparam.numneu + 1], p_neu[0])
            p_aneu = map(lambda x : x[sparam.numneu + 1], p_aneu[0])
            
            E_TeV = map(lambda x : x/param.TeV,E)
            
            plt.plot(E_TeV,p_neu,label = neu_label[ineu], lw = 4,linestyle = "dashed")
            plt.plot(E_TeV,p_aneu,label = aneu_label[ineu], lw = 4,linestyle = "dashed")
        
    
    plt.xlabel(r"$E_\nu \, [\mathrm{TeV}]$")
    plt.ylabel(r"$P_{absorption}$")
    
    #plt.loglog()
    plt.semilogx()
    
    plt.legend(loc = "upper right")
    
    plt.xlim(10.0e-3,1.0)
    plt.ylim(0.0,1.0)
    
    path = "../plots/"
    
    filename = "SunAbsorptionProbabilityMC_"+str(iineu)+".eps"
    
    plt.savefig(path+filename)
    
def PlotCompareSunAbsorptionProbabilityMCvsANA(iineu,Enumin,Enumax,param,datapath = "../data/SunOscProbabilities/"):
    """Plots absorption probability due to CC interactions using MC data.
    # Enumin        : minimum neutrino energy       [eV]
    # Enumax        : maximum neutrino energy       [eV]
    # param         : physics parameter
    """
    
    #E = gt.MidPoint(gt.LogSpaceEnergies(Enumin,Enumax,binnum = 200))
    neu_label  = ["$\\nu_e$","$\\nu_\\mu$","$\\nu_\\tau$"]
    aneu_label = ["$\\bar{\\nu}_e$","$\\bar{\\nu}_\\mu$","$\\bar{\\nu}_\\tau$"]
    
    ## setup canvas ##
    
    fig = plt.figure()
    ax = plt.subplot(111)
            
    mpl.rcParams['axes.labelsize'] = "large"
    mpl.rcParams['xtick.labelsize'] = "large"
    mpl.rcParams['legend.fontsize'] = "small"
            
    mpl.rcParams['font.size'] = 15
    
    ## end setup canvas ##
    
    ## FOR PARAM ##
     
    for ineu in [iineu]:
        
        filename_neutrino       = "DataNeuOscProb_RK_neutrino_Emin_"+str(Enumin/param.GeV)+"_GeV_Emax_"+str(Enumax/param.GeV)+"_GeV_ineu_"+str(ineu)+"_param_"+param.name+".dat"
        filename_antineutrino   = "DataNeuOscProb_RK_antineutrino_Emin_"+str(Enumin/param.GeV)+"_GeV_Emax_"+str(Enumax/param.GeV)+"_GeV_ineu_"+str(ineu)+"_param_"+param.name+".dat"
        
        file_neu  = open(datapath+filename_neutrino,'r')
        file_aneu = open(datapath+filename_antineutrino,'r')
        
        p_neu = []
        p_aneu = []   
        
        gt.hreadfilev4(file_neu,p_neu,param)
        gt.hreadfilev4(file_aneu,p_aneu,param)
        
        file_neu.close()
        file_aneu.close()
        
        E = map(lambda x : x[0], p_neu[0])
        p_neu  = map(lambda x : x[param.numneu + 1], p_neu[0])
        p_aneu = map(lambda x : x[param.numneu + 1], p_aneu[0])
        
        E_TeV = map(lambda x : x/param.TeV,E)
        
        plt.plot(E_TeV,p_neu,label = neu_label[ineu], lw = 4, linestyle = "solid")
        plt.plot(E_TeV,p_aneu,label = aneu_label[ineu], lw = 4, linestyle = "solid")
        
    ## BEGIN ANALITIC ##
    
    datapath = "../data/NeutrinoAbsorptionProb/"
      
    filename_neutrino       = "DatNeuSunAbsorptionProbability_1.0_GeV_1000000.0_GeV_neutrino.dat"
    filename_antineutrino   = "DatNeuSunAbsorptionProbability_1.0_GeV_1000000.0_GeV_antineutrino.dat"
    
    file_neu  = open(datapath+filename_neutrino,'r')
    file_aneu = open(datapath+filename_antineutrino,'r')
    
    p_neu = []
    p_aneu = []   
    
    gt.hreadfilev4(file_neu,p_neu,param)
    gt.hreadfilev4(file_aneu,p_aneu,param)
    
    file_neu.close()
    file_aneu.close()
    
    E = map(lambda x : x[0], p_neu[0])
    p_neu = map(lambda x : x[1], p_neu[0])
    p_aneu = map(lambda x : x[1], p_aneu[0])
    
    E_TeV = map(lambda x : x/param.TeV,E)
    
    plt.plot(E_TeV,p_neu,label = "DI "+neu_label[iineu], lw = 4)
    plt.plot(E_TeV,p_aneu,label = "DI "+aneu_label[iineu], lw = 4)
    
    ## END ANALITIC ##
    
    plt.xlabel(r"$E_\nu \, [\mathrm{TeV}]$")
    plt.ylabel(r"$P_{absorption}$")
    
    #plt.loglog()
    plt.semilogx()
    
    plt.legend(loc = "upper right")
    
    plt.xlim(10.0e-3,1.0)
    plt.ylim(0.0,1.0)
    
    path = "../plots/"
    
    filename = "PlotCompareSunAbsorptionProbabilityMC_"+str(iineu)+".eps"
    
    plt.savefig(path+filename)        

    

    
#===============================================================================
# Testing
#===============================================================================
    
if __name__ == '__main__':    
    Sun = bd.Sun()
    pc = PC.PhysicsConstants()
    
    pc.name = "STD"
    
    pcc = PC.PhysicsConstants()
    pcc.name = "3+2"
    
    Enumin = 10.0*pc.GeV
    Enumax = 10000.0*pc.GeV
    
    
    [PlotCompareSunAbsorptionProbabilityMCvsANA(iineu,Enumin,Enumax,pc,datapath = "../data/SunOscProbabilities/test/") for iineu in [0,1,2]]
    
    #[PlotSunAbsorptionProbabilityMC(iineu,Enumin,Enumax,pc,pcc,datapath = "../data/SunOscProbabilities/test/") for iineu in [0,1,2]]

    
    quit()
    #
    #pct = PC.PhysicsConstants()
    #pct.name = "3+3_v30b"
    #pct.numneu = 6
    #
    #thsterile = 0.087
    #dmsterile = 0.10
    #
    #pct.th14 = thsterile
    #pct.th25 = thsterile
    #pct.th36 = thsterile
    #
    #pct.dm41sq = dmsterile
    #pct.dm51sq = dmsterile
    #pct.dm61sq = dmsterile
    #
    #pct.Refresh()
    #
    #pc.name = "STD_v30b"
    #
    ##pct.name = "3+3"
    ##pc.name = "STD"
    #
    #PlotNewIcecubeXsectionLimit(pct,pc,None,True)
    #quit()
    #
    #
    #
    #
    ##GenerateConvolutedDMFlux(pc)
    ##quit()
    #
    #PlotIcecubeLimitCompareMCwithANA(pct)
    #quit()
    #
    ##datapath = "/home/carlos/workspace/Solarneu/data/ParameterScan/3+3v2/"
    ##
    ###CalculateParameterScan(datapath)
    ###quit()
    ##
    ##datapath = "/home/carlos/workspace/Solarneu/data/ParameterScan/3+3fv/"
    ##
    ##PlotParameterScan(datapath)
    ##quit()
    ##PlotOscProb(0,Enumin,Enumax,[pc,pct], datapath = "../data/SunOscProbabilities/test/")
    ##quit()
    
    
    pcc = PC.PhysicsConstants()
    pcc.name = "3+2"
    pcc.linestyle = "dashed"
    pcc.markerstyle = "o"
    pcc.colorstyle = "blue"
    pcc.numneu = 5 
    
    pcc.th12 = 0.601264
    pcc.th13 = 0.0
    pcc.th23 = 0.737177
        
    pcc.th14 = -0.148782
    pcc.th24 = 0.150234
    pcc.th34 = -0.0223698
    pcc.th15 = 0.123729
    pcc.th25 = 0.118667
    pcc.th34 = 0.15
    pcc.th35 = -0.00507371

    pcc.dm21sq = 7.6e-05
    pcc.dm31sq = 0.00238375    
    pcc.dm41sq = 0.900335
    pcc.dm51sq = 0.472633
    
    pcc.delta1 = 4.60779
    pcc.delta2 = 0.269928
    pcc.delta3 = 2.41591
    
    pcc.Refresh()
    
    pcc.neutype = "antineutrino"
    
    PlotOscProb(0,Enumin,Enumax,[pc,pcc], datapath = "../data/SunOscProbabilities/test/")
    quit()     
    
    #PlotSolarParameters()
    
    #quit()
    
    #pc.dm21sq = 8.1e-5
    #pc.dm32sq = 2.2e-3
    #pc.th12   = 33.0*np.pi/180.0
    #pc.th13   = 0.0
    #pc.th23   = np.pi/4.0
    #pc.Refresh()
    #PlotDMEventSpectraCompare('bb',50,Sun,pc)
    #PlotDMFluxSpectraCompare('bb',50,Sun,pc)
    #PlotDMEventSpectraMultiCh(50,Sun,pc)
    #PlotDMEventSpectraMultiCh(100,Sun,pc)
    #PlotDMFluxSpectraMultiCh(50,Sun,pc)
    #PlotDMFluxSpectraMultiCh(100,Sun,pc)
    
    #PlotDMFluxSpectraSingleCh('tautau',50,Sun,pc)
    #PlotDMFluxSpectraSingleCh('tautau',100,Sun,pc)
    #PlotDMFluxSpectraMultiCh(300,Sun,pc)
    #PlotDMEventSpectraMultiCh(50,Sun,pc)
    #PlotDMFluxSpectraMultiChCorrections(100,Sun,pc)
    
    #PlotDMMC_Icecube(pc)
    #PlotDMMC_MTonDetector(pc)
    #quit()
    
    #PlotDMFLuxProduction(100.0,'tautau')
    #quit()
    #PlotDMFLuxProduction(50.0,'tautau')
    
    #PlotDM_Icecube(pc)
    
    pcc = PC.PhysicsConstants()
    pcc.name = "3+2"
    pcc.linestyle = "dashed"
    pcc.markerstyle = "o"
    pcc.colorstyle = "blue"
    pcc.numneu = 5 
    
    pcc.th14 = 0.129599 # higher increases resonances amplitud # if <0.1 no effect
    pcc.th24 = 0.171848 # no effect on pee
    pcc.th15 = 0.138442
    pcc.th25 = 0.149991
    pcc.th34 = 0.15#0.05 # no effect on pee
    pcc.th13 = 0.0
    pcc.dm41sq = -0.47 # higher moves resonante to right / proportional
    pcc.dm51sq = -0.87
    pcc.delta1 = 0.0
    
    #pcc.neutype = "antineutrino"
    #pc.neutype = "antineutrino"
    
    pcc.Refresh()
    
    import hotshot as hs
    import hotshot.stats as hstat
    
    prof = hs.Profile("MC_profile_0.prof")
    
    datapath = "../data/myMC/"+pct.name+"/"+pct.neutype+"/"
    
    PlotNewIcecubeXsectionLimit(pct,pc,None,True)
    #quit()
    
    ch = 'bb'
    dm = 1.0*pct.TeV
    sig = 1.8e-41*pct.cm**2
    
    benchtime = prof.runcall(DM.DMFNeuFluxMCDetv2,ch,dm,sig,pct,False,datapath)
    
    #benchtime = prof.runcall(PlotNewIcecubeXsectionLimit,pct,pc,None,True)
    
    prof.close()

    stat = hstat.load("MC_profile_0.prof")
    stat.strip_dirs()
    stat.print_stats()        
    
    #quit()
    
    #pct = PC.PhysicsConstants()
    #pct.name = "3+3v2"
    #pct.numneu = 6
    #
    #thsterile = 0.087
    #dmsterile = 0.10
    #
    #pct.th14 = thsterile
    #pct.th25 = thsterile
    #pct.th36 = thsterile
    #
    #pct.dm41sq = dmsterile
    #pct.dm51sq = dmsterile
    #pct.dm61sq = dmsterile
    #
    #pct.Refresh()

    pc.neutype = "antineutrino"
    pcc.neutype = "antineutrino"    
    
    datapath = "../data/SunOscProbabilities/"
    
    PlotOscProb(2,Enumin,Enumax,[pc,pcc])
    quit()    
    
    Enumin = 1.0*pc.GeV
    Enumax = 1.0*pc.TeV
    
    E = 1.0*pc.TeV

    quit()
    PlotAnnihilationRatesLimit(pc)
    quit()
    PlotNewIcecubeLimit(pc)
    quit()
    PlotCompareProbabilitiesMC(pcc,pc)
    quit()
    #PlotFluxAtDetector(pc,pc)
    PlotDMMC_MTonDetector(pcc,pc)
    quit()
    
    #PlotOscProb(Enumin,Enumax,[pc,pcc])
    PlotSingleNeuCompositionCompare(E,Sun,pcc,pc)
    #PlotOscProb(Enumin,Enumax,[pc,pcc])
    
    quit()
    
    PMNS = no.mixmatrix(pc)
    fM2  = no.flavorM2(pc)
    
    Ri   = 0.05*Sun.Radius*pcc.km
    Rf   = Sun.Radius*pcc.km
#    #
#    ##PlotDMFluxSpectraMultiCh(100,Sun,pcc)
    E = 10*pcc.MeV#1012773328388.837
#    ##E    = 1.0*pcc.TeV#10.0*pcc.TeV
#    ##prob = [no.AdiabaticProbability(i,i,E,Ri,Sun,PMNS,fM2,pcc) for i in range(5)]
#    #
    #track = Sun.track(Ri,Ri,Rf)
#    #
#    #print no.AvgNeuProb_RK(1,E,track,Sun,fM2,pcc)
#    
#    def singlepoint():
#        no.AvgNeuProb_RK(1,E,track,Sun,fM2,pcc)
#    
#    
    import hotshot
    import hotshot.stats as hstat
    #prof = hotshot.Profile("RKprof.prof")
    #benchtime = prof.runcall(no.AvgNeuProb_RK(1,E,track,Sun,fM2,pcc))
    #prof.close()
    #
    #print benchtime
    #profile.runcall(bd.Sun())
    #quit()
#    profile.run('singlepoint()')
    #profile.run(no.AvgNeuProb_RK(1,E,track,Sun,fM2,pcc))
    
#    quit()
    ##print no.AvgNeuProb_RK(2,E,track,Sun,fM2,pcc)
    #
    #no.DataNeuOscProb(0,1.0*pcc.GeV,10*pcc.GeV,Sun,pcc)
    
    #benchtime = prof.runcall(no.DataNeuOscProb,0,1.0*pcc.GeV,10*pcc.GeV,Sun,pcc)
    #prof.close()
    
    stat = hstat.load("RKprof.prof")
    stat.strip_dirs()
    stat.print_stats()
    
    quit()
    #pc.Refresh()
    #print no.AdiabaticProbability(0,0,E,Ri,Sun,PMNS,fM2,pc),no.AdiabaticProbability(0,1,E,Ri,Sun,PMNS,fM2,pc),no.AdiabaticProbability(0,2,E,Ri,Sun,PMNS,fM2,pc)#,no.AdiabaticProbability(0,3,E,Ri,Sun,PMNS,fM2,pcc),no.AdiabaticProbability(0,4,E,Ri,Sun,PMNS,fM2,pcc)
    #print no.AdiabaticProbability(0,0,E,Ri,Sun,PMNS,fM2,pc)+no.AdiabaticProbability(0,1,E,Ri,Sun,PMNS,fM2,pc)+no.AdiabaticProbability(0,2,E,Ri,Sun,PMNS,fM2,pc)#+no.AdiabaticProbability(0,3,E,Ri,Sun,PMNS,fM2,pcc)+no.AdiabaticProbability(0,4,E,Ri,Sun,PMNS,fM2,pcc)
    #
    #quit()
    #print no.AdiabaticProbability(0,3,E,Ri,Sun,PMNS,fM2,pcc),no.AdiabaticProbability(0,4,E,Ri,Sun,PMNS,fM2,pcc)
    #print no.AdiabaticProbability(1,3,E,Ri,Sun,PMNS,fM2,pcc),no.AdiabaticProbability(1,4,E,Ri,Sun,PMNS,fM2,pcc)
    #print no.AdiabaticProbability(2,3,E,Ri,Sun,PMNS,fM2,pcc),no.AdiabaticProbability(2,4,E,Ri,Sun,PMNS,fM2,pcc)
    #
    #print prob
    #
    ##E = 10.0*pcc.GeV
    ##prob = [no.AdiabaticProbability(i,i,E,0.05,Sun,PMNS,fM2,pcc) for i in range(3)]
    ##print prob
    #
    #quit()
        
    #PlotNeuCompositionCompare(250*pc.GeV,Sun,pcc,pc)
    #quit()
    
    #PlotDM_Icecube(pcc,pc)
    PlotDM_MTonDetector(pcc,pc)
    quit()
    datapath = "../data/"
    
    import generaltools as gt
    
    file = open(datapath + "DataNeuOscProb_RK_antineutrino_Emin_1.0_GeV_Emax_1000.0_GeV_ineu_2_param_2+3.dat", 'r')
    h,d  = gt.hreadfilev2(file)
    
    E = np.array(d)[:,0]
    P = np.array(d)[:,1]
    
    plt.plot(E,P)
    
    plt.show()
    
    
