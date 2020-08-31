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
import neutrinocommon.neu.neuosc as no
import neutrinocommon.physconst.physicsconstants as PC
import neutrinocommon.tools.generaltools as gt
    
#===============================================================================
# STERILE NEUTRINO OSCILLATION PROBABILITY
#===============================================================================

def PlotNeuOscProb(ineu,fneu,Ri,Enumin,Enumax,body,param,rkcomp,plot_path = "../plots/"):
    """Plots P(neu_ineu -> neu_fneu) as a function of the Energy
    # ineu,fneu     : 0 (electron), 1 (muon), 2 (tau)
    # Ri            : neutrino production point in the body
    # Enumin        : minimum neutrino energy
    # Enumax        : maximum neutrino energy
    # body          : body where the neutrino propagates
    # param         : physics parameter set
    # rkcomp        : toggle RK validation
    # plot_path     : path to save plot             [string]
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

    filename = plot_path+"PlotNeuOscProb_Emin_"+str(Enumin/scale)+"_"+scalen+"_Emax_"+str(Enumax/scale)+"_"+scalen+"_nui_"+sineulabel+"_nuf_"+sfneulabel+".png"
    plt.savefig(filename)
    
def PlotNeuOscProbSun(ineu,Enumin,Enumax,param,rkcomp,plot_path = "../plots/"):
    """Plots P(neu_ineu -> neu_fneu) as a function of the Energy from an initial flavor state (ineu) to all final flavor states (fneu) on the sun
    # ineu,fneu     : 0 (electron), 1 (muon), 2 (tau)
    # Enumin        : minimum neutrino energy       [eV]
    # Enumax        : maximum neutrino energy       [eV]
    # param         : physics parameter set list    [param_1,param_2,...,param_n]
    # rkcomp        : toggle RK validation          [boolean]
    # plot_path     : path to save plot             [string]    
    # If rkcomp = True, will use RK to compute the probability and compare
    """
    body    = bd.Sun()
    
    fig = plt.figure(figsize = (4*3+2,6))
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
        
    neulabel    = {0 : "e",1 : "\\mu", 2 : "\\tau",3 : "{s_1}",4 : "{s_2}"}
    sneulabel   = {0 : "e",1 : "mu", 2 : "tau",3 : "s1",4 : "s2"} 
            
    # RK points
    ERKstep = (np.log10(Enumax)-np.log10(Enumin))/(20.0)    
    ERK = np.arange(np.log10(Enumin),np.log10(Enumax),ERKstep)
    ERK = map(lambda E : (10**E)/scale,ERK)
    
    Estep = (Enumax-Enumin)/npoints        
    Enu = np.arange(Enumin,Enumax,Estep)/scale
    # Generating plot
    for fneu in [0,1,2]:
        fig.add_subplot(1,3,fneu+1)
        try:
            for i,p in enumerate(param):
                p.Refresh()
                
                PMNS = no.mixmatrix(p)
                fM2 = no.flavorM2(p)
                Ri = 0.01*body.Radius*p.km
                Rf = body.Radius*p.km
                track = body.track(Ri,Ri,Rf)
                
                Pe = map(lambda E : no.AdiabaticProbability(ineu,fneu,E*scale,Ri,body,PMNS,fM2,p),Enu)
                
                plt.xlabel(r"$\mathrm{E}_\nu\mathrm{["+scalen+"]}$")
                plt.ylabel("$ P(\\nu_{"+neulabel[ineu]+"} \\rightarrow \\nu_{"+neulabel[fneu]+"}) $")
                plt.plot(Enu,Pe,color = 'black',linestyle = p.linestyle,label='$ P_{\\mathrm{Adb.}-\\mathrm{'+p.name+'}}(\\nu_'+neulabel[ineu]+'\\rightarrow \\nu_'+neulabel[fneu]+')$')
                
                if rkcomp :
                    PRK = map(lambda E : no.AvgNeuProb_RK_STD(ineu,fneu,E*scale,p),ERK)
                    plt.plot(ERK,PRK,p.markerstyle,color=p.colorstyle,markersize=6,markeredgecolor=p.colorstyle,markeredgewidth=1,label='$ P_{\\mathrm{RK}-\\mathrm{'+p.name+'}}(\\nu_'+neulabel[ineu]+'\\rightarrow \\nu_'+neulabel[fneu]+')$')
                        
        except (TypeError,AttributeError):
            param.Refresh()
            
            PMNS = no.mixmatrix(param)
            fM2 = no.flavorM2(param)
            Ri = 0.01*body.Radius*param.km
            Rf = body.Radius*param.km
            track = body.track(Ri,Ri,Rf)
            Pe = map(lambda E : no.AdiabaticProbability(ineu,fneu,E*scale,Ri,body,PMNS,fM2,param),Enu)
            
            plt.xlabel(r"$\mathrm{E}_\nu\mathrm{["+scalen+"]}$")
            plt.ylabel("$ P(\\nu_{"+neulabel[ineu]+"} \\rightarrow \\nu_{"+neulabel[fneu]+"}) $")
            plt.plot(Enu,Pe,color = 'black',linestyle = param.linestyle,label='$ P_{\\mathrm{Adb.}-\\mathrm{'+param.name+'}}(\\nu_'+neulabel[ineu]+'\\rightarrow \\nu_'+neulabel[fneu]+')$')
            
            if rkcomp :
                PRK = map(lambda E : no.probneu(ineu,fneu,E*scale,track,body,fM2,param),ERK)
                plt.plot(ERK,PRK,param.markerstyle,color=param.colorstyle,markersize=6,markeredgecolor=param.colorstyle,markeredgewidth=1,label='$ P_{\\mathrm{RK}-\\mathrm{'+param.name+'}}(\\nu_'+neulabel[ineu]+'\\rightarrow \\nu_'+neulabel[fneu]+')$')
        
        plt.semilogx()
        plt.axis([Enumin/scale,Enumax/scale,0.0,1.0])
    
    try:
        plt.suptitle(param[0].neutype+" oscillations probabilities")
    except TypeError:
        plt.suptitle(param.neutype+" oscillations probabilities")
    fig.subplots_adjust(left=0.05, right=0.85,wspace = 0.35, top = 0.85, bottom = 0.15)
    mpl.rcParams['legend.fontsize'] = "small"
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,fancybox = True)
#    path = "../plots/"
    try:
        filename = plot_path+"PlotNeuOscProbSun_Emin_"+str(Enumin/scale)+"_"+scalen+"_Emax_"+str(Enumax/scale)+"_"+scalen+"_nui_"+neulabel[ineu]+"_"+param[0].neutype+".png"
    except TypeError:
        filename = plot_path+"PlotNeuOscProbSun_Emin_"+str(Enumin/scale)+"_"+scalen+"_Emax_"+str(Enumax/scale)+"_"+scalen+"_nui_"+neulabel[ineu]+"_"+param.neutype+".png"
    plt.savefig(filename)    

def PlotEigenvaluesRho(E,body,param):
    plt.figure()
    param.Refresh()
    fM2 = no.flavorM2(param)
    if(E/param.MeV <= 500.0) :
        scale =  param.MeV
        scalen = "MeV"
    elif(E/param.GeV <= 500.0) :
        scale =  param.GeV
        scalen = "GeV"
    else :
        scale =  param.TeV
        scalen = "TeV"
    R = np.arange(1.0,0.01,-0.001)
    Rho = map(lambda r : body.rdensity(r), R)
    EE = map(lambda x : no.Eeigenvals(E, x, body, fM2, param),R)
    plt.xlabel(r"$\rho \mathrm{[g/cm^{-3}]}$")
    plt.ylabel(r"$\mathrm{E[eV]}$")
    plt.plot(Rho,EE,color = 'black')
    #plt.loglog()
    plt.semilogx()
    #plt.ylim(1.0e-14,1.0e-12)

    filename = "PlotEnergyEigenvaluesRho_E_"+str(E/scale)+"_"+scalen+".png"
    plt.savefig(plot_path + filename)
    
#===============================================================================
# NEUTRINO EIGENVALUES PLOT
#===============================================================================     
    
def PlotEigenvaluesEnergy(Emin,Emax,xRi,body,param,plot_path = "../plots/"):
    plt.figure()
    param.Refresh()
    fM2 = no.flavorM2(param)
    ordmag = np.log10(Emax)-np.log10(Emin)
    npoints = 1000.0*ordmag
    # begin energy scale # 
    if(Emax/param.MeV <= 500.0) :
            scale =  param.MeV
            scalen = "[MeV]"
    elif(Emax/param.GeV <= 500.0) :
            scale =  param.GeV
            scalen = "[GeV]"
    else :
            scale =  param.TeV
            scalen = "[TeV]"
    # end energy scale #
    Estep = (Emax-Emin)/npoints
    Enu = np.arange(Emin,Emax,Estep)/scale
    EE = map(lambda E : no.Eeigenvals(E*scale, xRi, body, fM2, param),Enu)
    plt.xlabel(r"$E_\nu \mathrm{"+scalen+"}$")
    plt.ylabel(r"$\mathrm{E[eV]}$")
    plt.plot(Enu,EE,color = 'black')
    #plt.loglog()
    plt.semilogx()
    plt.savefig(plot_path+"PlotEnergyEigenvaluesEnergy.png")
    
def PlotEigenvaluesRadius(E,body,param,plot_path = "../plots/"):
    plt.figure()
    R = np.arange(0.01,1.0,0.01)
    param.Refresh()
    fM2 = no.flavorM2(param)
    EE = map(lambda x : no.Eeigenvals(E, x, body, fM2, param),R)
    plt.xlabel(r"$\frac{R}{R_\odot}$")
    plt.ylabel(r"$\mathrm{E[eV]}$")
    plt.plot(R,EE,color = 'black')
    plt.semilogy()
    plt.savefig(plot_path+"PlotEnergyEigenvaluesRadius.png")
    
#===============================================================================
# NEUTRINO STATES COMPOSITION
#===============================================================================    

def PlotNeuComposition(E,body,param,plot_path = "../plots/",xlim = None, fmt = "png"):
    """ Plots the composition of neutrinos as a function of mass states, flavors, and density.
    
    E        :    neutrino energy [eV]
    body     :    body with the asociated density profile.
    param    :    set of physical parameters used to make the plot. param can be a list. In this case both sets will be plotted/label
    """
    fig = plt.figure(figsize = (4*param.numneu,8))
    #B Initializing variables
    param.Refresh()
    fM2 = no.flavorM2(param)
    R = np.arange(1.0,0.01,-0.001)
    Rho = map(lambda r : body.rdensity(r), R)
    #E Initializing variables
    #B Estimating Energy Scale
    if(E/param.MeV <= 500.0) :
        scale =  param.MeV
        scalen = "MeV"
    elif(E/param.GeV <= 500.0) :
        scale =  param.GeV
        scalen = "GeV"
    else :
        scale =  param.TeV
        scalen = "TeV"
    #E Estimating Energy Scale
    #B Adding title
    
    tit = "Energy : "+str(E/scale)+" "+scalen+ " Parameters :"#+" $\\th_{12}$ = " + str(param.th12) + " $\\th_{23}$ = " + str(param.th23) + " $\\th_{13}$ = "+str(param.th13)
    atit = []
    try:
        [[ atit.append(" $\\theta_{"+str(j)+str(i)+"}$ = "+format(param.th[j][i],'.4f')) for i in range(1,param.numneu+1) if i>j] for j in range(1,param.numneu+1) ]
        [ atit.append(" $\\Delta m^2_{"+str(j)+str(1)+"}$ = "+format(param.dm2[1][j],'.3e')) for j in range(2,param.numneu+1) ]
    except TypeError:
        [[ atit.append(" $\\theta_{"+str(j)+str(i)+"}$ = "+"{:.4}".format(param.th[j][i])) for i in range(1,param.numneu+1) if i>j] for j in range(1,param.numneu+1) ]
        [ atit.append(" $\\Delta m^2_{"+str(j)+str(1)+"}$ = "+"{:.3}".format(param.dm2[1][j])) for j in range(2,param.numneu+1) ]
        
    for i in range(len(atit)):
        tit = tit + atit[i]
    plt.suptitle(tit,horizontalalignment='center')
    #E Adding title
    ##B PLOTTING MASS BASIS AS FUNCTION OF FLAVOR BASIS
    for i in np.arange(0,param.numneu,1):
        fig.add_subplot(2,param.numneu,i+1)
        flavor = False
        NeuComp = map(lambda x : no.NeuComposition(i,E, x, body, fM2, param,flavor),R)
        plt.xlabel(r"$\rho \mathrm{[g/cm^{3}]}$")
        pp = plt.plot(Rho,NeuComp,lw = '3')#,label =("$v_e$","$v_\mu$"))
        #Solar density
        ps = plt.vlines(150, 0.0, 1.0, linestyle = "dashed", label = r"$\rho_S$")
        #B plt format
        plt.title(r"Composition of $\nu_"+str(i+1)+"$")
        plt.semilogx()
        plt.ylim(0.0,1.0)
        plt.yticks(np.arange(0.0,1.1,0.1))
        
        if xlim == None :
            plt.xlim(Rho[0],Rho[-1])      
        else : 
            plt.xlim(xlim)           
        #plt.xscale()
        if i == param.numneu - 1 :
            plots = [] 
            for e in pp :
                plots.append(e)
            plots.append(ps) 
            leg = ["$\\nu_e$","$\\nu_\mu$","$\\nu_\\tau$"]
            ss =  ["$\\nu_{s"+str(i)+"}$" for i in np.arange(1,param.numneu-3+1,1)]
            if ss != []:
                leg.extend(ss)
            leg.append("$\\rho_\odot$")
            plt.legend(plots,leg,bbox_to_anchor = (1.5, 0.90))
        fig.subplots_adjust(left=0.05, right=0.85,hspace = 0.35)#,wspace = 0.25, top = 0.95, bottom = 0.05)
        #E plt format
    ##E PLOTTING MASS BASIS AS FUNCTION OF FLAVOR BASIS
    ##B PLOTTING FLAVOR BASIS AS FUNCTION OF MASS BASIS    
    for i in np.arange(param.numneu,2*param.numneu,1):
        fig.add_subplot(2,param.numneu,i+1)
        flavor = True
        NeuComp = map(lambda x : no.NeuComposition(i-param.numneu,E, x, body, fM2, param,flavor),R)
        plt.xlabel(r"$\rho \mathrm{[g/cm^{3}]}$")
        pp = plt.plot(Rho,NeuComp,lw = '3')#,label =("$v_e$","$v_\mu$"))
        #Solar density
        ps = plt.vlines(150, 0.0, 1.0, linestyle = "dashed", label = r"$\rho_S$")
        #B plt format
        if i == param.numneu :
            plt.title(r"Composition of $\nu_e$")
        elif i == param.numneu + 1:
            plt.title(r"Composition of $\nu_\mu$")
        elif i == param.numneu + 2:
            plt.title(r"Composition of $\nu_\tau$")
        else : 
            plt.title(r"Composition of $\nu_{s"+str(i+1-param.numneu-3)+"}$")
            
        plt.semilogx()
        plt.ylim(0.0,1.0)
        
        if xlim == None :
            plt.xlim(Rho[0],Rho[-1])      
        else : 
            plt.xlim(xlim)                
#        plt.xlim(0.0,1000.0)
        plt.yticks(np.arange(0.0,1.1,0.1))
        #plt.xscale()
        if i == 2*param.numneu - 1 :
            plots = [] 
            for e in pp :
                plots.append(e)
            plots.append(ps) 
            leg =  ["$\\nu_{"+str(i)+"}$" for i in np.arange(1,param.numneu+1,1)]
            leg.append("$\\rho_\odot$")
            plt.legend(plots,leg,bbox_to_anchor = (1.5, 0.90))
        fig.subplots_adjust(left=0.05, right=0.85)#, top = 0.85, bottom = 0.05)        
    ##E PLOTTING FLAVOR BASIS AS FUNCTION OF MASS BASIS
#    path = "../plots/"
    filename = "PlotNeuComposition_E_"+str(E/scale)+"_"+scalen+"_"+param.name+"."+fmt
    
    plt.savefig(plot_path + filename)
    

def PlotNeuCompositionCompare(E,body,param,sparam = PC.PhysicsConstants(),plot_path = "../plots/",xlim = None, fmt = "png"):
    """ Plots the composition of neutrinos as a function of mass states, flavors, and density. Compares with STD.
    
    E        :    neutrino energy [eV]
    body     :    body with the asociated density profile.
    param    :    set of physical parameters used to make the plot. param can be a list.
    sparam   :    standard parameters
    """
    fig = plt.figure(figsize = (4*param.numneu,8))
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
    #B Initializing variables
    param.Refresh()
    fM2 = no.flavorM2(param)
    fM2STD = no.flavorM2(sparam)
    R = np.arange(1.0,0.01,-0.001)
    Rho = map(lambda r : body.rdensity(r), R)
    #E Initializing variables
    #B Estimating Energy Scale
    if(E/param.MeV <= 500.0) :
        scale =  param.MeV
        scalen = "MeV"
    elif(E/param.GeV <= 500.0) :
        scale =  param.GeV
        scalen = "GeV"
    else :
        scale =  param.TeV
        scalen = "TeV"
    #E Estimating Energy Scale
    #B Adding title
    
    tit = "Energy : "+str(E/scale)+" "+scalen+ " Parameters :"#+" $\\th_{12}$ = " + str(param.th12) + " $\\th_{23}$ = " + str(param.th23) + " $\\th_{13}$ = "+str(param.th13)
    atit = []
    try:
        [[ atit.append(" $\\theta_{"+str(j)+str(i)+"}$ = "+format(param.th[j][i],'.4f')) for i in range(1,param.numneu+1) if i>j] for j in range(1,param.numneu+1) ]
        [ atit.append(" $\\Delta m^2_{"+str(j)+str(1)+"}$ = "+format(param.dm2[1][j],'.3e')) for j in range(2,param.numneu+1) ]
    except TypeError:
        [[ atit.append(" $\\theta_{"+str(j)+str(i)+"}$ = "+"{:.4}".format(param.th[j][i])) for i in range(1,param.numneu+1) if i>j] for j in range(1,param.numneu+1) ]
        [ atit.append(" $\\Delta m^2_{"+str(j)+str(1)+"}$ = "+"{:.3}".format(param.dm2[1][j])) for j in range(2,param.numneu+1) ]
        
    for i in range(len(atit)):
        tit = tit + atit[i]
    plt.suptitle(tit,horizontalalignment='center')
    
    #E Adding title
    ##B PLOTTING MASS BASIS AS FUNCTION OF FLAVOR BASIS
    for i in np.arange(0,param.numneu,1):
        fig.add_subplot(2,param.numneu,i+1)
        flavor = False
        NeuComp = map(lambda x : no.NeuComposition(i,E, x, body, fM2, param,flavor),R)
        plt.xlabel(r"$\rho \mathrm{[g/cm^{3}]}$")
        pp = []
        for k in range(param.numneu):
            kNeuComp = map(lambda x: x[k],NeuComp)
            pp.append(plt.plot(Rho,kNeuComp, color = colors[k],lw = 3))
        if i<=2 :
            NeuCompSTD = map(lambda x : no.NeuComposition(i,E, x, body, fM2STD, sparam,flavor),R)
            for k in range(sparam.numneu):
                kNeuCompSTD = map(lambda x: x[k],NeuCompSTD)
                plt.plot(Rho,kNeuCompSTD, color = colors[k], linestyle = 'dashed',lw = 3)
        #Solar density
        ps = plt.vlines(150, 0.0, 1.0, linestyle = "dashed", label = r"$\rho_S$")
        #B plt format
        plt.title(r"Composition of $\nu_"+str(i+1)+"$")
        plt.semilogx()
        plt.ylim(0.0,1.0)
        plt.xlim(0.0,1000.0)
        plt.yticks(np.arange(0.0,1.1,0.1))
        #plt.xscale()
        if xlim == None :
            plt.xlim(Rho[0],Rho[-1])      
        else : 
            plt.xlim(xlim)            
        if i == param.numneu - 1 :
            plots = [] 
            for e in pp :
                plots.append(e[0])
            plots.append(ps)
            leg = ["$\\nu_e$","$\\nu_\mu$","$\\nu_\\tau$"]
            ss =  ["$\\nu_{s"+str(i)+"}$" for i in np.arange(1,param.numneu-3+1,1)]
            if ss != []:
                leg.extend(ss)
            leg.append("$\\rho_\odot$")
            plt.legend(plots,leg,bbox_to_anchor = (1.5, 0.90))
        fig.subplots_adjust(left=0.05, right=0.85,hspace = 0.35)#,wspace = 0.25, top = 0.95, bottom = 0.05)
        #E plt format
    ##E PLOTTING MASS BASIS AS FUNCTION OF FLAVOR BASIS
    ##B PLOTTING FLAVOR BASIS AS FUNCTION OF MASS BASIS    
    for i in np.arange(param.numneu,2*param.numneu,1):
        fig.add_subplot(2,param.numneu,i+1)
        flavor = True
        NeuComp = map(lambda x : no.NeuComposition(i-param.numneu,E, x, body, fM2, param,flavor),R)
        plt.xlabel(r"$\rho \mathrm{[g/cm^{-3}]}$")
        pp = []
        for k in range(param.numneu):
            kNeuComp = map(lambda x: x[k],NeuComp)
            pp.append(plt.plot(Rho,kNeuComp, color = colors[k],lw = 3))
        if i<=2+param.numneu :
            NeuCompSTD = map(lambda x : no.NeuComposition(i-param.numneu,E, x, body, fM2STD, sparam,flavor),R)
            for k in range(sparam.numneu):
                kNeuCompSTD = map(lambda x: x[k],NeuCompSTD)
                plt.plot(Rho,kNeuCompSTD, color = colors[k], linestyle = 'dashed',lw = 3)
        #Solar density
        ps = plt.vlines(150, 0.0, 1.0, linestyle = "dashed", label = r"$\rho_S$")
        #B plt format
        if i == param.numneu :
            plt.title(r"Composition of $\nu_e$")
        elif i == param.numneu + 1:
            plt.title(r"Composition of $\nu_\mu$")
        elif i == param.numneu + 2:
            plt.title(r"Composition of $\nu_\tau$")
        else : 
            plt.title(r"Composition of $\nu_{s"+str(i+1-param.numneu-3)+"}$")
            
        plt.semilogx()
        plt.ylim(0.0,1.0)
        plt.yticks(np.arange(0.0,1.1,0.1))
        
        if xlim == None :
            plt.xlim(Rho[0],Rho[-1])      
        else : 
            plt.xlim(xlim)    
        #plt.xscale()
        if i == 2*param.numneu - 1 :
            plots = [] 
            for e in pp :
                plots.append(e[0])
            plots.append(ps) 
            leg =  ["$\\nu_{"+str(i)+"}$" for i in np.arange(1,param.numneu+1,1)]
            leg.append("$\\rho_\odot$")
            plt.legend(plots,leg,bbox_to_anchor = (1.5, 0.90))
        fig.subplots_adjust(left=0.05, right=0.85)#, top = 0.85, bottom = 0.05)        
    ##E PLOTTING FLAVOR BASIS AS FUNCTION OF MASS BASIS
    plt.suptitle("*Dashed colored lines are 3-flavor standard oscillations.", x = 0.15, y = 0.03)

    filename = "PlotNeuComposition_E_"+str(E/scale)+"_"+scalen+"_"+param.name+"."+fmt
    plt.savefig(plot_path + filename)
    
#===============================================================================
# Compare Prob Osc
#===============================================================================    

def PlotCompareProbabilitiesMC(param,sparam = PC.PhysicsConstants(),plot_path = "../plots/"):
    """ Plot osc. probabilities from MC compared with analitic calculation of probabilities.
    
    @type   param   :   physicsconstants
    @param  param   :   set of physical parameters used to make the plot.
    @type   sparam  :   physicsconstants
    @param  sparam  :   standard parameters
    
    @rtype          :   plot
    @return         :   generates probability comparison plot
    """
    fig = plt.figure()
    
    Enumin  = 1.0
    Enumax  = 1.0e3
    Enustep = 10.0
    
    Enu = np.arange(Enumin,Enumax,Enustep)

    # param1

    param.name  =  "2+3"
    param.neutype = "neutrino"
    datapath = "../data/myMC/"+param.name+"/"+param.neutype+"/"

    inter_nu, inter_anu = DM.DMOscProbabilitiesMC(param,False, datapath = datapath, crosscheck = False)    
    mc_osc_prob  = [float(inter_nu(E)) for E in Enu]
    
    # sparam
    sparam.name = "STD"
    param.neutype = "neutrino"
    sdatapath = "../data/myMC/"+sparam.name+"/"+sparam.neutype+"/"
    
    inter_nu, inter_anu = DM.DMOscProbabilitiesMC(sparam,False, datapath = sdatapath, crosscheck = False)    
    smc_osc_prob  = [float(inter_nu(E)) for E in Enu]
    
    # analitic
    
    ana_osc_prob = [no.AvgNeuProb_RK_STD(0,1,E*param.GeV,param)+no.AvgNeuProb_RK_STD(1,1,E*param.GeV,param)+no.AvgNeuProb_RK_STD(2,1,E*param.GeV,param) for E in Enu]
    
    plt.plot(Enu,mc_osc_prob,label = 'MC-2+3')
    plt.plot(Enu,smc_osc_prob,label = 'MC-STD')
    plt.plot(Enu,ana_osc_prob,label = 'ANA')
    
    plt.legend()
    
#    path = "../plots/"
    filename = "PlotCompareProbabilitiesMC.png"
    plt.savefig(plot_path+filename)                
    

#===============================================================================
# Adiabaticity
#===============================================================================

def PlotAdiabaticityParameter(param,sparam = PC.PhysicsConstants()):
    """ Plots the adiabaticity parameter from for a given energy range on the sun.
    
    @type   param   :   physicsconstants
    @param  param   :   set of physical parameters used to make the plot.
    @type   sparam  :   physicsconstants
    @param  sparam  :   standard parameters
    
    @rtype          :   plot
    @return         :   generates the spectra plot
    """
    
    import Scientific.Functions.Derivatives as sfd
    
    
    Sun = bd.Sun()
    
    print Sun.rdensity(sfd.DerivVar(0.5))
    
    ## th_24 resonance search
    
    E = np.arange(1.0,1000.0,10.0)
    
    adiabaticity_s1 = (np.sin(2.0*param.th24)**2/np.cos(2.0*param.th24))*(param.dmsq41/(2.0*E*param.GeV))*()
    adiabaticity_s2 = (np.sin(2.0*param.th24)**2/np.cos(2.0*param.th24))*(param.dmsq51/(2.0*E*param.GeV))*()

        
#===============================================================================
# Probability Running time
#===============================================================================    
    
def PlotRKProbabilityTime(param):
    """ Plots the CPU time in second to calculate the probability.
    
    @type   param   :   physicsconstants
    @param  param   :   set of physical parameters used to make the plot.
    
    @rtype          :   plot
    @return         :   generates the spectra plot
    """
    import hotshot as hs
    import hotshot.stats as hstat
    import re
    prof = hs.Profile("RKprob.prof")
    
    E = np.arange(1.0,1000.0,50.0)
    
    Sun = bd.Sun()
    Ri   = 0.01*Sun.Radius*param.km
    Rf   = Sun.Radius*param.km
    track = Sun.track(Ri,Ri,Rf)
    fM2  = no.flavorM2(param)
    
    EE = 1.0*param.GeV
    
    benchtime = prof.runcall(no.AvgNeuProb_RK,1,EE,track,Sun,fM2,param)
    prof.close()
    print benchtime
    
    stat = hstat.load("RKprob.prof")
    stat.strip_dirs()
    
    stat.dump_stats('RKprob.stat')
    
    #re.match('calls',stat.print_stats(0))
    stat.print_stats()       

#===============================================================================
# PUBLICATION FINAL PLOTS
#=============================================================================== 

def PlotOscProb(iineu,Enumin,Enumax,param, datapath = "../data/SunOscProbabilities/",plotpath = "../plots",plot_survival_probability = False,filename_append = '', fmt = 'eps'):
    """Plots P(neu_ineu -> neu_fneu) as a function of the Energy from an initial flavor state (ineu)
    to all final flavor states (fneu) on the sun
    # iineu         : 0 (electron), 1 (muon), 2 (tau)
    # Enumin        : minimum neutrino energy       [eV]
    # Enumax        : maximum neutrino energy       [eV]
    # param         : physics parameter set list    [param_1,param_2,...,param_n]
    """
    plt.cla()
    plt.clf()
    plt.close()    
    #fig = plt.figure(figsize = (4*3+2,6))
    fig = plt.figure(figsize = (10,7.5))
    ax = plt.subplot(111)
    
    mpl.rcParams['axes.labelsize'] = "xx-large"
    mpl.rcParams['xtick.labelsize'] = "xx-large"
    mpl.rcParams['legend.fontsize'] = "small"
    
    mpl.rcParams['font.size'] = 30
    
    ordmag = np.log10(Enumax)-np.log10(Enumin)
    npoints = 1000.0*ordmag
    # Figuring out best energy scale
    try:
        if(Enumax/param[0].MeV <= 500.0) :
            scale =  param[0].MeV
            scalen = "MeV"
        elif(Enumax/param[0].GeV <= 1000.0) :
            scale =  param[0].GeV
            scalen = "GeV"
        else :
            scale =  param[0].GeV#param[0].TeV
            scalen = "GeV"#"TeV"
    except (TypeError,AttributeError):
        if(Enumax/param.MeV <= 500.0) :
            scale =  param.MeV
            scalen = "MeV"
        elif(Enumax/param.GeV <= 1000.0) :
            scale =  param.GeV
            scalen = "GeV"
        else :
            scale =  param.TeV
            scalen = "TeV"
            
    try : 
        Emin = Enumin/param[0].GeV
        Emax = Enumax/param[0].GeV
    except :
        Emin = Enumin/param.GeV
        Emax = Enumax/param.GeV        
        
    
    neulabel    = {0 : "e",1 : "\\mu", 2 : "\\tau",3 : "{s_1}",4 : "{s_2}",5 : "{s_3}"}
    sneulabel   = {0 : "e",1 : "mu", 2 : "tau",3 : "s1",4 : "s2",5 : "s3"} 
            
    # RK points
    #ERKstep = (np.log10(Enumax)-np.log10(Enumin))/(20.0)    
    #ERK = np.arange(np.log10(Enumin),np.log10(Enumax),ERKstep)
    #ERK = map(lambda E : (10**E)/scale,ERK)
    #ERK.append(1000.0)
    
    ERK = gt.MidPoint(gt.LogSpaceEnergies(Enumin/scale,Enumax/scale, binnum = 200))
    
    ERK = [ERK[i] for i in range(len(ERK)-1)]
    
    ERK = [ERK[i] for i in range(len(ERK)-1)]
    
    Estep = (Enumax-Enumin)/npoints        
    Enu = np.arange(Enumin,Enumax,Estep)/scale
    # Generating plot
    
    #totalPRK = [0.0]*len(ERK)
    
    #colors   = ['orange', 'r', 'k', 'c', 'm', 'y', 'k']
    colors   = ['b', 'r', 'g','c', 'm', 'y', 'k']
    linestyles = ['--','-.',':','-..','-']
    
    for ineu,fneu in [[iineu,0],[iineu,1],[iineu,2],[iineu,5]]:#[[0,0],[1,1],[2,2],[2,1],[2,0],[0,0]]:
        for i,p in enumerate(param):
            p.Refresh()
            plt.xlabel(r"$\mathrm{E}_\nu\mathrm{["+scalen+"]}$")
            plt.ylabel("$ \mathrm{Probability}$")
            
            if (p.name == "STD" or p.name == "STD_XXX" or p.numneu == 3) and fneu > 2:
                pass
            else :
                if fneu == 5:
                    fneu = 3
                    PRK_3 = map(lambda E : float(no.InterOscProb(ineu,fneu,E*scale,p,datapath,Emin,Emax,filename_append = filename_append)),ERK)
                    fneu = 4
                    PRK_4 = map(lambda E : float(no.InterOscProb(ineu,fneu,E*scale,p,datapath,Emin,Emax,filename_append = filename_append)),ERK)
                    if p.numneu > 5 :
                        fneu = 5
                        PRK_5 = map(lambda E : float(no.InterOscProb(ineu,fneu,E*scale,p,datapath,Emin,Emax,filename_append = filename_append)),ERK)
                    else :
                        PRK_5 = map(lambda E : 0.0*E ,ERK)
                    #PRK = [PRK_3[i] + PRK_4[i] for i in range(len(PRK_3))]
                    PRK = [PRK_3[i] + PRK_4[i] + PRK_5[i] for i in range(len(PRK_3))]
                    
                    if p.neutype == "neutrino":
                        plt.plot(ERK,PRK,linestyle = linestyles[-1],label='$ \mathrm{P}(\\nu_'+neulabel[ineu]+'\\rightarrow \\nu_s)$',color = p.colorstyle, lw = 6,solid_joinstyle = 'bevel')
                    elif p.neutype == "antineutrino":
                        plt.plot(ERK,PRK,linestyle = linestyles[-1],label='$ \mathrm{P}(\\bar{\\nu}_'+neulabel[ineu]+'\\rightarrow \\bar{\\nu}_s)$',color = p.colorstyle, lw = 6,solid_joinstyle = 'bevel')
                    else:
                        print "Wrong neutrino type."
                        quit()
                else :
                    PRK = map(lambda E : no.InterOscProb(ineu,fneu,E*scale,p,datapath,Emin,Emax,filename_append = filename_append),ERK)
                    if p.name == "STD" or p.name == "STD_XXX" or p.numneu == 3:
                        plt.plot(ERK,PRK,linestyle = linestyles[fneu],color = p.colorstyle, lw = 4)
                    else :
                        if p.neutype == "neutrino":
                            plt.plot(ERK,PRK,linestyle = linestyles[fneu],label='$ \mathrm{P}(\\nu_'+neulabel[ineu]+'\\rightarrow \\nu_'+neulabel[fneu]+')$', color = p.colorstyle, lw = 4)
                        elif p.neutype == "antineutrino":
                            plt.plot(ERK,PRK,linestyle = linestyles[fneu],label='$ \mathrm{P}(\\bar{\\nu}_'+neulabel[ineu]+'\\rightarrow \\bar{\\nu}_'+neulabel[fneu]+')$', color = p.colorstyle, lw = 4)
                        else:
                            print "Wrong neutrino type."
                            quit()                                
                
    if plot_survival_probability :
        for p in param:
            P_surival = map(lambda E : no.InterOscProb(iineu,p.numneu,E*scale,p,datapath,Emin,Emax,filename_append = filename_append),ERK)
            plt.plot(ERK,P_surival,linestyle = 'solid', lw = 4, color = p.colorstyle,solid_joinstyle = 'bevel')
        
    plt.semilogx()
    #plt.loglog()
    
    plt.axis([Enumin/scale,Enumax/scale,0.0,1.0])
    
    #fig.subplots_adjust(left=0.05, right=0.8,wspace = 0.35, top = 0.85, bottom = 0.15)
    #mpl.rcParams['legend.fontsize'] = "small"
    #plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,fancybox = True)
    plt.legend(loc='upper right',fancybox = True)
    
    ######################### begin extra legend ##########################
    for i,p in enumerate(param):
        # removing everything after an _ (underscore)
        try:
            paramname = re.search('.*(?=_)',p.name).group(0)
        except AttributeError:
            paramname = p.name
        
        # create text box
        boxt = osb.TextArea(paramname, textprops=dict(color="k"))
        boxd = osb.DrawingArea(60, 20, 0, 0)
        el = ptc.Ellipse((10, 10), width=5, height=5, angle=0, fc=p.colorstyle, edgecolor = 'none')
        boxd.add_artist(el)
        
        box = osb.HPacker(children=[boxt, boxd],
              align="center",
              pad=0, sep=1)
        
        # anchor boxes
        anchored_box = osb.AnchoredOffsetbox(loc=2,
                            child=box, pad=0.25,
                            frameon=False,
                            bbox_to_anchor=(0.0, 1.0-0.06*i),
                            bbox_transform=ax.transAxes,
                            borderpad=0.,
                            )
        
        ax.add_artist(anchored_box)
    ########################## end extra legend ##############################
    
    fig.subplots_adjust(bottom = 0.12, top = 0.95, left = 0.12, right = 0.95)
    
    path = plotpath
    
    mpl.rcParams['font.size'] = 30
    
    try:
        filename = path+"PlotOscProbability_ineu_"+str(iineu)
        for p in param:
            filename = filename + +"_" + p.name+"_"+p.neutype
        filename = filename + "."+fmt
    except TypeError:
        filename = path+"PlotOscProbability_ineu_"+str(iineu)+"_"+param.name+"_"+param.neutype+"."+fmt
        
    plt.savefig(filename, dpi = 1200)
    
    plt.clf()
    plt.close()

def PlotSingleNeuCompositionCompare(E,body,param,sparam = PC.PhysicsConstants()):
    """ Plots the composition of a single mass neutrino state.
    
    E        :    neutrino energy [eV]
    body     :    body with the asociated density profile.
    param    :    set of physical parameters used to make the plot. param can be a list.
    sparam   :    standard parameters
    """
    fig = plt.figure()
    ax = plt.subplot(111)
    
    mpl.rcParams['axes.labelsize'] = "x-large"
    mpl.rcParams['xtick.labelsize'] = "x-large"
    mpl.rcParams['legend.fontsize'] = "small"
    
    mpl.rcParams['font.size'] = 18    
    
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
    linestyles = ['--','-.',':','-','-']
    
    #B Initializing variables
    param.Refresh()
    fM2 = no.flavorM2(param)
    sparam.Refresh()
    fM2STD = no.flavorM2(sparam)
    R = np.arange(1.0,0.01,-0.001)
    Rho = map(lambda r : body.rdensity(r), R)
    #E Initializing variables
    #B Estimating Energy Scale
    if(E/param.MeV <= 500.0) :
        scale =  param.MeV
        scalen = "MeV"
    elif(E/param.GeV <= 500.0) :
        scale =  param.GeV
        scalen = "GeV"
    else :
        scale =  param.TeV
        scalen = "TeV"
    #E Estimating Energy Scale
    
    #B Adding title
    #tit = "Energy : "+str(E/scale)+" "+scalen+ " Parameters :"#+" $\\th_{12}$ = " + str(param.th12) + " $\\th_{23}$ = " + str(param.th23) + " $\\th_{13}$ = "+str(param.th13)
    #atit = []
    #[[ atit.append(" $\\theta_{"+str(j)+str(i)+"}$ = "+format(param.th[j][i],'.4f')) for i in range(1,param.numneu+1) if i>j] for j in range(1,param.numneu+1) ]
    #[[ atit.append(" $\\Delta m^2_{"+str(i)+str(j)+"}$ = "+format(param.dm2[j][i],'.4f')) for i in range(1,param.numneu+1) if i>j and j == 1] for j in range(1,param.numneu+1) ]
    ##[[ atit.append(" $\\Delta m^2_{"+str(j)+str(i)+"}$ = "+format(param.dm2[j][i],'.4f')) for i in range(1,param.numneu+1) if i>j and j == 1] for j in range(1,param.numneu+1) ]
    #for i in range(len(atit)):
    #    tit = tit + atit[i]
    #plt.suptitle(tit,horizontalalignment='center')    
    #E Adding title
    
    ##B PLOTTING MASS BASIS AS FUNCTION OF FLAVOR BASIS
    for i in [1]:
        #fig.add_subplot(2,param.numneu,i+1)
        flavor = False
        NeuComp = map(lambda x : no.NeuComposition(i,E, x, body, fM2, param,flavor),R)
        plt.xlabel(r"$\rho \mathrm{[g/cm^{3}]}$")
        pp = []
        for k in range(param.numneu):
            kNeuComp = map(lambda x: x[k],NeuComp)
            
            # log interpolator
            rholog = gt.LogSpaceEnergies(float(Rho[0]),float(Rho[-1]),100)
            
            #print rholog , float(Rho[0]),float(Rho[-1])
            rholog[-1] = Rho[-1]
            
            inter_neu = interpolate.interp1d(Rho,kNeuComp)
            logkNeuComp = map(inter_neu,rholog)
            
            if k == 3:
                #pp.append(plt.plot(Rho,kNeuComp,'o-',color = 'r',markevery = 10,markeredgewidth = 0.0, ms = 2.0))
                pp.append(plt.plot(rholog,logkNeuComp,'x-',color = 'r',markevery = 10,markeredgewidth = 0.0, ms = 2.0,aa = True,solid_joinstyle = 'bevel'))
            elif k == 4:
                pp.append(plt.plot(Rho,kNeuComp,'o-',color = 'r',markevery = 10,markeredgewidth = 0.0, ms = 2.0))
            else :
                pp.append(plt.plot(Rho,kNeuComp,linestyle = linestyles[k] ,color = 'r'))
        if i<=2 :
            NeuCompSTD = map(lambda x : no.NeuComposition(i,E, x, body, fM2STD, sparam,flavor),R)
            for k in range(sparam.numneu):
                kNeuCompSTD = map(lambda x: x[k],NeuCompSTD)
                plt.plot(Rho,kNeuCompSTD, color = 'k', linestyle = linestyles[k])
        #Solar density
        #ps = plt.vlines(150, 0.0, 1.0, linestyle = "dashed", label = r"$\rho_S$")
        #B plt format
        plt.title(r"Composition of $\nu_"+str(i+1)+"$")
        plt.semilogx()
        plt.ylim(0.0,1.0)
        plt.xlim(1.0,150.0)
        plt.yticks(np.arange(0.0,1.1,0.1))
        xtt = [1.0,5.0,10.0,30.0,100.0]#,150.0]
        #plt.xticks(xtt)
        ax.set_xticks(xtt)
	ax.set_xticklabels(['$1$','$5$','$10$','$30$','$100$'])#,'$\\rho_\\odot = 150$'])
        
        #plt.xscale()
        #B LEGEND
        plots = [] 
        for e in pp :
            plots.append(e[0])
        #plots.append(ps)
        leg = ["$\\nu_e$","$\\nu_\mu$","$\\nu_\\tau$"]
        ss =  ["$\\nu_{s"+str(i)+"}$" for i in np.arange(1,param.numneu-3+1,1)]
        if ss != []:
            leg.extend(ss)
        leg = plt.legend(plots,leg,loc = 6,fancybox=True,bbox_to_anchor = (0.05, 0.75))
        leg.get_frame().set_alpha(0.25)
        #E LEGEND            
        #E plt format
        
        #B EXTRA LEGEND
        box1t = osb.TextArea("STD", textprops=dict(color="k"))
        box1d = osb.DrawingArea(60, 20, 0, 0)
        el1 = ptc.Ellipse((10, 10), width=5, height=5, angle=0, fc="k", edgecolor = 'none')
        box1d.add_artist(el1)
        
        box2t = osb.TextArea("2+3", textprops=dict(color="k"))
        box2d = osb.DrawingArea(60, 20, 0, 0)
        el2 = ptc.Ellipse((10, 10), width=5, height=5, angle=0, fc="r", edgecolor = 'none')
        box2d.add_artist(el2)
        
        box1 = osb.HPacker(children=[box1t, box1d],
                  align="center",
                  pad=0, sep=1)
        
        box2 = osb.HPacker(children=[box2t, box2d],
                  align="center",
                  pad=0, sep=5)
        
        anchored_box1 = osb.AnchoredOffsetbox(loc=9,
                                     child=box1, pad=5.0,
                                     frameon=False,
                                     #bbox_to_anchor=(0., 1.02),
                                     #bbox_transform=ax.transAxes,
                                     borderpad=0.,
                                     )
        
        anchored_box2 = osb.AnchoredOffsetbox(loc=9,
                                     child=box2, pad=6.0,
                                     frameon=False,
                                     #bbox_to_anchor=(0., 1.02),
                                     #bbox_transform=ax.transAxes,
                                     borderpad=0.,
                                     )    
        
        ax.add_artist(anchored_box1)
        ax.add_artist(anchored_box2)
        #E EXTRA LEGEND
        
    ##E PLOTTING MASS BASIS AS FUNCTION OF FLAVOR BASIS

    #plt.suptitle("*Dashed colored lines are 3-flavor standard oscillations.", x = 0.15, y = 0.03)
    path = "../plots/"
    filename = "PlotNeuComposition_E_"+str(E/scale)+"_"+scalen+"_FIG2.eps"
    plt.savefig(path + filename, dpi = 1200)    
    
#===============================================================================
# Testing
#===============================================================================
    
if __name__ == '__main__':    
	pass
