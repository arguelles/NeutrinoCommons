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
# DM EVENT SPECTRA AT DETECTOR
#===============================================================================        
    
    
def PlotDMEventSpectraCompare(ch,DMm,body,param,sparam = PC.PhysicsConstants()):
    """ Plot DM event spectra after propagation at Earth.
    
    Plots the spectra of the DM at production of a given I{annihilation channel} [ch] and
    dark matter mass. Also calculates and plots the spectra of the produced neutrinos
    at detection under standard 3-generation neutrino oscillation and non standard
    sterile neutrino oscillation. Propagation of the neutrino is made under the B{adiabatic hipotesis}.
    
    @type   ch      :   string
    @param  ch      :   dark matter annhilation channel
    @type   DMm     :   float
    @param  DMm     :   dark matter mass
    @type   body    :   body
    @param  body    :   body with the asociated density profile.
    @type   param   :   physicsconstants
    @param  param   :   set of physical parameters used to make the plot.
    @type   sparam  :   physicsconstants
    @param  sparam  :   standard parameters
    
    @rtype          :   plot
    @return         :   generates the spectra plot
    """

    fig = plt.figure(figsize = (4*3+1.5,4.2))
    
    Enu = np.arange(1.0,DMm,DMm/1000.0)
    Emu = map(lambda EE: (1.0-xs.ymuCC(EE,0))*EE,Enu)

    for flavor in [0,1,2]:
        fig.add_subplot(1,3,flavor+1)
        p = DM.DMParameters(flavor)
        DMspectraEarth = [DM.DMneuDet(flavor,E*param.GeV,ch,DMm*param.GeV,body,pc,True)*(E*np.log(10.0))/DMm for E in Enu]
        DMspectraEarth_NoOsc = [DM.DMneuDet(flavor,E*param.GeV,ch,DMm*param.GeV,body,pc,False)*(E*np.log(10.0))/DMm for E in Enu]
        EnuODMm = map(lambda x: x/DMm,Enu)
    
        plt.plot(Emu,DMspectraEarth, label = r"Std. 3-gen Osc.", linestyle = "solid", color = "black")
        plt.plot(Emu,DMspectraEarth_NoOsc, label = r"No Osc.", linestyle = "dashed", color = "black")
        
        plt.xlabel(r"$E_\nu$")
        plt.ylabel(r"$\frac{dN}{d\mathrm{log}_{10}(E_\nu)}$")
        plt.semilogx()
        plt.xlim(xmin = 3.0, xmax = DMm)
        plt.ylim(ymin = 0.0)
        if flavor == 0 :
            plt.title(r"$\nu_e$")
        elif flavor == 1:
            plt.title(r"$\nu_\mu$")
        elif flavor == 2:
            plt.title(r"$\nu_\tau$")
    fig.subplots_adjust(left=0.07, right=0.85,wspace = 0.3, top = 0.85, bottom = 0.15)
    mpl.rcParams['legend.fontsize'] = "small"
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,fancybox = True)
    plt.suptitle(r"$\nu$ from $\chi \chi \rightarrow " + DM.ChannelLabel(ch) + "$ with $\chi_{mass} = "+str(DMm)+" \mathrm{GeV}$")
    path = "../plots/"
    filename = "PlotDMEventSpectraChannel_"+ch+"_DMm_"+str(DMm)+"_GeV.png"
    plt.savefig(path + filename)
    
def PlotDMEventSpectraMultiCh(DMm,body,param,sparam = PC.PhysicsConstants()):
    """ Plot DM event spectra after propagation at Earth for all channels.
    
    Plots the spectra of the DM at production of a given I{annihilation channel} [ch] and
    dark matter mass. Also calculates and plots the spectra of the produced neutrinos
    at detection under standard 3-generation neutrino oscillation and non standard
    sterile neutrino oscillation. Propagation of the neutrino is made under the B{adiabatic hipotesis}.
    
    @type   ch      :   string
    @param  ch      :   dark matter annhilation channel
    @type   DMm     :   float
    @param  DMm     :   dark matter mass
    @type   body    :   body
    @param  body    :   body with the asociated density profile.
    @type   param   :   physicsconstants
    @param  param   :   set of physical parameters used to make the plot.
    @type   sparam  :   physicsconstants
    @param  sparam  :   standard parameters
    
    @rtype          :   plot
    @return         :   generates the spectra plot
    """

    fig = plt.figure(figsize = (4*3+1.5,4.2))
    channels = {'bb' : 0, 'tautau' : 1}
    #channels = {'bb' : 0, 'tautau' : 1, 'cc': 2}
    chname = ['$b\\bar{b}$','$\\tau\\bar{\\tau}$']
    #chname = ['$b\\bar{b}$','$\\tau\\bar{\\tau}$','$c\\bar{c}$']
    colors = ['orange', 'r', 'k', 'c', 'm', 'y', 'k']
    
    Enu = np.arange(1.0,DMm,DMm/1000.0)
    Emu = map(lambda EE: (1.0-xs.ymuCC_avg(EE,0))*EE,Enu)

    for flavor in [0,1,2]:
        fig.add_subplot(1,3,flavor+1)
        for ch,chid in channels.iteritems():
            DMspectraEarth          = [DM.DMneuDet(flavor,E*param.GeV,ch,DMm*param.GeV,body,param,True)*(E*np.log(10.0))/DMm*(1.0-xs.ymuCC_avg(E,0)) for E in Enu]
            DMspectraEarth_NoOsc    = [DM.DMneuDet(flavor,E*param.GeV,ch,DMm*param.GeV,body,param,False)*(E*np.log(10.0))/DMm*(1.0-xs.ymuCC_avg(E,0)) for E in Enu]
        
            plt.plot(Emu,DMspectraEarth, label = chname[chid]+" std. 3-gen Osc.", linestyle = "solid", color = colors[chid])
            plt.plot(Emu,DMspectraEarth_NoOsc, label = chname[chid]+" no Osc.", linestyle = "dashed", color = colors[chid])
        
        plt.xlabel(r"$E_\mu$")
        plt.ylabel(r"$\frac{dN}{d\mathrm{log}_{10}(E_\mu)}$")
        plt.semilogx()
        plt.xlim(xmin = 3.0, xmax = DMm)
        plt.ylim(ymin = 0.0)
        if flavor == 0 :
            plt.title(r"$\nu_e$")
        elif flavor == 1:
            plt.title(r"$\nu_\mu$")
        elif flavor == 2:
            plt.title(r"$\nu_\tau$")
    fig.subplots_adjust(left=0.07, right=0.85,wspace = 0.3, top = 0.85, bottom = 0.15)
    mpl.rcParams['legend.fontsize'] = "small"
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,fancybox = True)
    plt.suptitle(r"$\chi_{mass} = "+str(DMm)+" \mathrm{GeV}$")
    path = "../plots/"
    filename = "PlotDMEventSpectraMultiChannel_DMm_"+str(DMm)+"_GeV.png"
    plt.savefig(path + filename)    

#===============================================================================
# DM FLUX SPECTRA AT DETECTOR
#===============================================================================          
    
def PlotDMFluxSpectraCompare(ch,DMm,body,param,sparam = PC.PhysicsConstants()):
    """ Plot DM flux spectra after propagation at Earth.
    
    Plots the spectra of the DM at production of a given I{annihilation channel} [ch] and
    dark matter mass. Also calculates and plots the spectra of the produced neutrinos
    at detection under standard 3-generation neutrino oscillation and non standard
    sterile neutrino oscillation. Propagation of the neutrino is made under the B{adiabatic hipotesis}.
    
    @type   ch      :   string
    @param  ch      :   dark matter annhilation channel
    @type   DMm     :   float
    @param  DMm     :   dark matter mass
    @type   body    :   body
    @param  body    :   body with the asociated density profile.
    @type   param   :   physicsconstants
    @param  param   :   set of physical parameters used to make the plot.
    @type   sparam  :   physicsconstants
    @param  sparam  :   standard parameters
    
    @rtype          :   plot
    @return         :   generates the spectra plot
    """

    fig = plt.figure(figsize = (4*3+1.5,4.2))
    
    Enu = np.arange(1.0,DMm,DMm/1000.0)
    for flavor in [0,1,2]:
        fig.add_subplot(1,3,flavor+1)
        
        DMspectraEarth = [(DM.DMFluxneuDet(flavor,E*param.GeV,ch,DMm*param.GeV,body,pc,True)*(E*np.log(10.0))/DMm)*(param.meter**2*param.sec) for E in Enu]
        DMspectraEarth_NoOsc = [(DM.DMFluxneuDet(flavor,E*param.GeV,ch,DMm*param.GeV,body,pc,False)*(E*np.log(10.0))/DMm)*(param.meter**2*param.sec) for E in Enu]       
        plt.plot(Enu,DMspectraEarth, label = r"Std. 3-gen Osc.", linestyle = "solid", color = "black")
        plt.plot(Enu,DMspectraEarth_NoOsc, label = r"No Osc.", linestyle = "dashed", color = "black")
        
        plt.xlabel(r"$E_\nu$")
        plt.ylabel(r"$\frac{d\Phi}{d\mathrm{log}_{10}(E_\nu)} / \mathrm{m}^2 \mathrm{s}$")
        plt.semilogx()
        plt.xlim(xmin = 3.0, xmax = DMm)
        plt.ylim(ymin = 0.0)
        if flavor == 0 :
            plt.title(r"$\nu_e$")
        elif flavor == 1:
            plt.title(r"$\nu_\mu$")
        elif flavor == 2:
            plt.title(r"$\nu_\tau$")
    fig.subplots_adjust(left=0.07, right=0.80,wspace = 0.3, top = 0.85, bottom = 0.15)
    mpl.rcParams['legend.fontsize'] = "small"
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,fancybox = True)
    plt.suptitle(r"$\nu$ from $\chi \chi \rightarrow " + DM.ChannelLabel(ch) + "$ with $\chi_{mass} = "+str(DMm)+" \mathrm{GeV}$")
    path = "../plots/"
    filename = "PlotDMFluxSpectraChannel_"+ch+"_DMm_"+str(DMm)+"_GeV.png"
    plt.savefig(path + filename)

def PlotDMFluxSpectraSingleCh(ch,DMm,body,param,sparam = PC.PhysicsConstants()):
    """ Plot DM flux spectra after propagation at Earth in all channels.
    
    Plots the spectra of the DM at production of a given I{annihilation channel} [ch] and
    dark matter mass. Also calculates and plots the spectra of the produced neutrinos
    at detection under standard 3-generation neutrino oscillation and non standard
    sterile neutrino oscillation. Propagation of the neutrino is made under the B{adiabatic hipotesis}.
    
    @type   DMm     :   float
    @param  DMm     :   dark matter mass
    @type   body    :   body
    @param  body    :   body with the asociated density profile.
    @type   param   :   physicsconstants
    @param  param   :   set of physical parameters used to make the plot.
    @type   sparam  :   physicsconstants
    @param  sparam  :   standard parameters
    
    @rtype          :   plot
    @return         :   generates the spectra plot
    """

    #fig = plt.figure(figsize = (8,4.2))
    fig = plt.figure()
    colors = ['orange', 'r', 'k', 'c', 'm', 'y', 'k']
    
    Enu = np.arange(1.0,DMm,DMm/1000.0)
    Z   = map(lambda x: x/DMm, Enu)

    for flavor in [0,1,2]:
            
        DMspectraEarth_NoOsc = [(DM.DMFluxneuDet(flavor,E*param.GeV,ch,DMm*param.GeV,body,pc,False)*(E*np.log(10.0))/DMm) for E in Enu]
            
        param.neutype = "neutrino"
        DMspectraEarth_neutrino         = [(DM.DMFluxneuDet(flavor,E*param.GeV,ch,DMm*param.GeV,body,pc,True))*(param.cm**2)*1.0e35 for E in Enu]
        param.neutype = "antineutrino"
        DMspectraEarth_antineutrino     = [(DM.DMFluxneuDet(flavor,E*param.GeV,ch,DMm*param.GeV,body,pc,True))*(param.cm**2)*1.0e35 for E in Enu]
            
        # combine neutrino-antineutrino spectra
        DMspectra_combined = [(2.0*DMspectraEarth_neutrino[i] + DMspectraEarth_antineutrino[i])/3.0 for i in range(len(DMspectraEarth_neutrino))]
        
        if flavor == 0 :
            plt.plot(Z,DMspectra_combined, label = r"$\nu_e$", linestyle = "solid", color = colors[flavor])
        elif flavor == 1:
            plt.plot(Z,DMspectra_combined, label = r"$\nu_\mu$", linestyle = "solid", color = colors[flavor])
        elif flavor == 2:
            plt.plot(Z,DMspectra_combined, label = r"$\nu_\tau$", linestyle = "solid", color = colors[flavor])
        
        #plt.plot(Enu,DMspectraEarth_NoOsc, label = chname[chid]+" no Osc.", linestyle = "dashed", color = colors[chid])
        
        plt.xlabel(r"$z =E_\nu / \chi_m$")
        #plt.ylabel(r"$\frac{d\Phi}{d\mathrm{log}_{10}(E_\nu)} [\mathrm{m}^{-2} \mathrm{s}^{-1}]$")
        plt.ylabel(r"$\frac{d\Phi}{d z} [10^{-35}\mathrm{cm}^{-2}\mathrm{ann}^-1]$")
        #plt.semilogx()
        plt.semilogy()
        
        plt.xlim(xmin = 0.0, xmax = 1.0)
        plt.ylim(ymin = 0.0)

    fig.subplots_adjust(left=0.07, right=0.85, wspace = 0.3, top = 0.85, bottom = 0.15)
    mpl.rcParams['legend.fontsize'] = "small"
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,fancybox = True)
    plt.suptitle(r"$\nu$ from $\chi \chi \rightarrow " + DM.ChannelLabel(ch) + "$ $\chi_{mass} = "+str(DMm)+" \mathrm{GeV}$")
    #plt.suptitle(r"$\Phi = (2\Phi_\nu + \Phi_{\bar{\nu}})/3$",x = 0.925, y = 0.35)
    path = "../plots/"
    filename = "PlotDMFluxSpectraSingleChannel_"+ch+"_DMm_"+str(DMm)+"_GeV_perann.png"
    plt.savefig(path + filename)
    
    
def PlotDMFluxSpectraMultiCh(DMm,body,param,sparam = PC.PhysicsConstants()):
    """ Plot DM flux spectra after propagation at Earth in all channels.
    
    Plots the spectra of the DM at production of a given I{annihilation channel} [ch] and
    dark matter mass. Also calculates and plots the spectra of the produced neutrinos
    at detection under standard 3-generation neutrino oscillation and non standard
    sterile neutrino oscillation. Propagation of the neutrino is made under the B{adiabatic hipotesis}.
    
    @type   DMm     :   float
    @param  DMm     :   dark matter mass
    @type   body    :   body
    @param  body    :   body with the asociated density profile.
    @type   param   :   physicsconstants
    @param  param   :   set of physical parameters used to make the plot.
    @type   sparam  :   physicsconstants
    @param  sparam  :   standard parameters
    
    @rtype          :   plot
    @return         :   generates the spectra plot
    """

    fig = plt.figure(figsize = (4*3+1.5,4.2))
    #channels = DM.ch
    channels = {'bb' : 0, 'tautau' : 1, 'cc': 2}
    chname = ['$b\\bar{b}$','$\\tau\\bar{\\tau}$','$c\\bar{c}$']
    colors = ['orange', 'r', 'k', 'c', 'm', 'y', 'k']
    
    Enu = np.arange(1.0,DMm,DMm/1000.0)

    for flavor in [0,1,2]:
        fig.add_subplot(1,3,flavor+1)
        
        for ch,chid in channels.iteritems():
            #DMspectraEarth_NoOsc = [(DM.DMFluxneuDet(flavor,E*param.GeV,ch,DMm*param.GeV,body,pc,False)*(E*np.log(10.0))/DMm)*(param.meter**2*param.sec) for E in Enu]
            
            #param.neutype = "neutrino"
            #DMspectraEarth_neutrino         = [(DM.DMFluxneuDet(flavor,E*param.GeV,ch,DMm*param.GeV,body,pc,True)*(E*np.log(10.0))/DMm)*(param.meter**2*param.sec) for E in Enu]
            #param.neutype = "antineutrino"
            #DMspectraEarth_antineutrino     = [(DM.DMFluxneuDet(flavor,E*param.GeV,ch,DMm*param.GeV,body,pc,True)*(E*np.log(10.0))/DMm)*(param.meter**2*param.sec) for E in Enu]
            
            DMspectraEarth_NoOsc = [(DM.DMFluxneuDet(flavor,E*param.GeV,ch,DMm*param.GeV,body,pc,False)*(E*np.log(10.0))/DMm) for E in Enu]
            
            param.neutype = "neutrino"
            DMspectraEarth_neutrino         = [(DM.DMFluxneuDet(flavor,E*param.GeV,ch,DMm*param.GeV,body,pc,True))*(param.cm**2)*1.0e35 for E in Enu]
            param.neutype = "antineutrino"
            DMspectraEarth_antineutrino     = [(DM.DMFluxneuDet(flavor,E*param.GeV,ch,DMm*param.GeV,body,pc,True))*(param.cm**2)*1.0e35 for E in Enu]
            
            # combine neutrino-antineutrino spectra
            DMspectra_combined = [(2.0*DMspectraEarth_neutrino[i] + DMspectraEarth_antineutrino[i])/3.0 for i in range(len(DMspectraEarth_neutrino))]
            
            plt.plot(Enu,DMspectra_combined, label = chname[chid]+" std. 3-gen Osc.", linestyle = "solid", color = colors[chid])
            #plt.plot(Enu,DMspectraEarth_NoOsc, label = chname[chid]+" no Osc.", linestyle = "dashed", color = colors[chid])
        
        plt.xlabel(r"$E_\nu$")
        #plt.ylabel(r"$\frac{d\Phi}{d\mathrm{log}_{10}(E_\nu)} [\mathrm{m}^{-2} \mathrm{s}^{-1}]$")
        plt.ylabel(r"$\frac{d\Phi}{d z} [10^{-35}\mathrm{cm}^{-2}\mathrm{ann}^-1]$")
        #plt.semilogx()
        plt.semilogy()
        
        plt.xlim(xmin = 3.0, xmax = DMm)
        plt.ylim(ymin = 0.0)
        if flavor == 0 :
            plt.title(r"$\nu_e$")
        elif flavor == 1:
            plt.title(r"$\nu_\mu$")
        elif flavor == 2:
            plt.title(r"$\nu_\tau$")
    fig.subplots_adjust(left=0.07, right=0.85, wspace = 0.3, top = 0.85, bottom = 0.15)
    mpl.rcParams['legend.fontsize'] = "small"
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,fancybox = True)
    plt.suptitle("$\chi_{mass} = "+str(DMm)+" \mathrm{GeV}$")
    plt.suptitle(r"$\Phi = (2\Phi_\nu + \Phi_{\bar{\nu}})/3$",x = 0.925, y = 0.35)
    path = "../plots/"
    filename = "PlotDMFluxSpectraMultiChannel_DMm_"+str(DMm)+"_GeV_perann.png"
    plt.savefig(path + filename)
    
def PlotDMFluxSpectraMultiChCorrections(DMm,body,param,sparam = PC.PhysicsConstants()):
    """ Plot DM flux spectra after propagation at Earth comparing no osc. with osc.
    
    Plots the spectra of the DM at production of a given I{annihilation channel} [ch] and
    dark matter mass. Also calculates and plots the spectra of the produced neutrinos
    at detection under standard 3-generation neutrino oscillation and non standard
    sterile neutrino oscillation. Propagation of the neutrino is made under the B{adiabatic hipotesis}.
    
    @type   DMm     :   float
    @param  DMm     :   dark matter mass
    @type   body    :   body
    @param  body    :   body with the asociated density profile.
    @type   param   :   physicsconstants
    @param  param   :   set of physical parameters used to make the plot.
    @type   sparam  :   physicsconstants
    @param  sparam  :   standard parameters
    
    @rtype          :   plot
    @return         :   generates the spectra plot
    """

    fig = plt.figure(figsize = (4*3+1.5,4.2))
    #channels = DM.ch
    channels = {'bb' : 0, 'tautau' : 1, 'cc': 2}
    chname = ['$b\\bar{b}$','$\\tau\\bar{\\tau}$','$c\\bar{c}$']
    colors = ['orange', 'r', 'k', 'c', 'm', 'y', 'k']
    
    Enu = np.arange(1.0,DMm,DMm/1000.0)
    for flavor in [0,1,2]:
        fig.add_subplot(1,3,flavor+1)
        
        for ch,chid in channels.iteritems():
            DMspectraEarth = [(DM.DMFluxneuDet(flavor,E*param.GeV,ch,DMm*param.GeV,body,pc,True)*(E*np.log(10.0))/DMm)*(param.meter**2*param.sec) for E in Enu]
            DMspectraEarth_NoOsc = [(DM.DMFluxneuDet(flavor,E*param.GeV,ch,DMm*param.GeV,body,pc,False)*(E*np.log(10.0))/DMm)*(param.meter**2*param.sec) for E in Enu]
            spectracorrection = [DMspectraEarth[i]/DMspectraEarth_NoOsc[i] for i in range(len(DMspectraEarth_NoOsc))]
            plt.plot(Enu,spectracorrection, label = chname[chid], linestyle = "solid", color = colors[chid])
        
        plt.xlabel(r"$E_\nu \mathrm{[GeV]}$")
        plt.ylabel("Correction due to propagation")
        plt.loglog()       
        plt.xlim(xmin = 3.0, xmax = DMm)
        plt.ylim(ymin = 0.0,ymax = 10.0)
        plt.yticks([0.1,0.3,1.0,3.0,10.0])
        if flavor == 0 :
            plt.title(r"$\nu_e$")
        elif flavor == 1:
            plt.title(r"$\nu_\mu$")
        elif flavor == 2:
            plt.title(r"$\nu_\tau$")
    fig.subplots_adjust(left=0.07, right=0.85,wspace = 0.3, top = 0.85, bottom = 0.15)
    mpl.rcParams['legend.fontsize'] = "small"
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,fancybox = True)
    plt.suptitle("$\chi_{mass} = "+str(DMm)+" \mathrm{GeV}$")
    path = "../plots/"
    filename = "PlotDMFluxSpectraMultiChannelCorrection_DMm_"+str(DMm)+"_GeV.png"
    plt.savefig(path + filename)            

#===============================================================================
# DM FLUX SPECTRA AT PRODUCTION
#===============================================================================          
    
def PlotDMFLuxProduction(DMm,ch):
    """ Plots DM pseudoflux for a given channel.
    @type  DMm          :      float
    @param DMm          :      dark matter mass [GeV]
    @type  DMm          :      float
    @param DMm          :      dark matter mass [GeV]
    
    @rtype              :      plot
    @return             :      generates plot
    """
    fig = plt.figure(figsize = (6,4))
    pc = PC.PhysicsConstants()
    
    Enu = np.arange(0.0,DMm,DMm/100.0)
    #p = DM.DMParameters(0)
    #DMPFlux_e = map(lambda E : DM.DMFlux(E*pc.GeV,DMm*pc.GeV,ch,p)*(E*np.log(10.0)/DMm),Enu)
    #p = DM.DMParameters(2)
    #DMPFlux_t = map(lambda E : DM.DMFlux(E*pc.GeV,DMm*pc.GeV,ch,p)*(E*np.log(10.0)/DMm),Enu)
    
    #p = DM.DMParameters(0)
    #DMPFlux_e = map(lambda E : DM.DMFlux(E*pc.GeV,DMm*pc.GeV,ch,p),Enu)
    #p = DM.DMParameters(2)
    #DMPFlux_t = map(lambda E : DM.DMFlux(E*pc.GeV,DMm*pc.GeV,ch,p),Enu)
    
    DMPFlux_e = map(lambda E : DM.DMSweFlux(E,0,ch,DMm),Enu)
    DMPFlux_m = map(lambda E : DM.DMSweFlux(E,2,ch,DMm),Enu)
    DMPFlux_t = map(lambda E : DM.DMSweFlux(E,4,ch,DMm),Enu)
    
    plt.plot(Enu,DMPFlux_e, label = r"$\nu_e$ at production")
    plt.plot(Enu,DMPFlux_m, label = r"$\nu_\mu$ at production")
    plt.plot(Enu,DMPFlux_t, label = r"$\nu_\tau$ at production")
    plt.title("Neutrino event number at production")
    plt.ylabel(r"$\frac{\mathrm{dN}}{\mathrm{dz}} [\mathrm{ann}^{-1}]$")
    #plt.ylabel(r"$\frac{\mathrm{dN}}{d\mathrm{log_{10}}E_\nu}$")
    plt.xlabel(r"$E_\nu $")
    #plt.xlabel(r"$E_\nu / \chi_m$")
    #plt.semilogx()
    plt.semilogy()
    mpl.rcParams['legend.fontsize'] = "small"
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,fancybox = True)
    fig.subplots_adjust(left=0.12, right=0.65,wspace = 0.3, top = 0.85, bottom = 0.15)
    path = "../plots/"
    filename = "PlotDMFluxProduction_DMm_"+str(DMm)+"_GeV_ch_"+ch+".png"
    plt.savefig(path+filename)
    
#===============================================================================
# MUON EVENT SPECTRA AT ICECUBE-LIKE DETECTOR
#===============================================================================


def PlotDM_Icecube(param,sparam = PC.PhysicsConstants()):
    """ Plot muon event spectra from DM annihilation at the Sun on selected
    channels. Can compare to non-standard neutrino oscillation scenarios.
    Uses Icecube simulation on L{ice.NevtSpectra}. Considering a 1 year
    exposure time.
    
    @type   param   :   physicsconstants
    @param  param   :   set of physical parameters used to make the plot.
    @type   sparam  :   physicsconstants
    @param  sparam  :   standard parameters
    
    @rtype          :   plot
    @return         :   generates the spectra plot
    """
    years       = 1.0
    
    DM_mass     = [50.0*sparam.GeV,100.0*sparam.GeV,1.0*sparam.TeV]  # DARK MATTER MASS
    channels    = {'bb' : 0, 'tautau' : 1}#, 'cc': 2}
    chname      = ['$b\\bar{b}$','$\\tau\\bar{\\tau}$','$c\\bar{c}$']
    colors      = ['orange', 'r', 'k', 'c', 'm', 'y', 'k']
    
    fig = plt.figure(figsize = (3*len(DM_mass)+2,4))
    
    for i,DMm in enumerate(DM_mass):
        fig.add_subplot(1,len(DM_mass),i+1)
        for ch,chid in channels.iteritems():
            
            ## STANRDARD PARAMETER SET
            sparam.Refresh()
            
            # define muon neutrinos flux from DM annhiliation
            DM_flux_STD         = DM.DMFluxAtDetector(ch,DMm,sparam)
            Enu                 = np.arange(DM_flux_STD.minenergy/param.GeV,DM_flux_STD.maxenergy/param.GeV,0.1)
            
            # calculate muon event spectra at icecube-like detector
            if sparam.neutype == "neutrino":
                evt_spectra_STD = ice.NevtSpectra(DM_flux_STD,0,years,sparam)
                Emu             = map(lambda EE : ice.MuonEnergy(EE*sparam.GeV,0,sparam)/sparam.GeV, Enu)
            elif sparam.neutype == "antineutrino":
                evt_spectra_STD = ice.NevtSpectra(DM_flux_STD,1,years,sparam)
                Emu             = map(lambda EE : ice.MuonEnergy(EE*sparam.GeV,1,sparam)/sparam.GeV, Enu)
            else :
                print "Wrong neutrino type."
                quit()
                
            # converting to "italian" scale
            #evt_spectra_logmuE_STD = []
            evt_spectra_muE_STD = []
            for i,E in enumerate(Enu):
                #evt_spectra_logmuE_STD.append(np.log(10.0)*E*evt_spectra_STD[i])
                evt_spectra_muE_STD.append((E/Emu[i])*evt_spectra_STD[i])
                                        
            # plot
            #plt.plot(Emu,evt_spectra_logmuE_STD,label = chname[chid], linestyle = "solid", color = colors[chid])
            plt.plot(Emu,evt_spectra_muE_STD,label = chname[chid]+" "+sparam.name, linestyle = "solid", color = colors[chid])
            
            ## OTHER PARAMETER SET
            param.Refresh()
            
            # define muon neutrinos flux from DM annhiliation
            DM_flux_STD         = DM.DMFluxAtDetector(ch,DMm,param)
            Enu                 = np.arange(DM_flux_STD.minenergy/sparam.GeV,DM_flux_STD.maxenergy/sparam.GeV,0.1)
            
            # calculate muon event spectra at icecube-like detector
            if sparam.neutype == "neutrino":
                evt_spectra_STD = ice.NevtSpectra(DM_flux_STD,0,years,sparam)
                Emu             = map(lambda EE : ice.MuonEnergy(EE*sparam.GeV,0,sparam)/sparam.GeV, Enu)
            elif sparam.neutype == "antineutrino":
                evt_spectra_STD = ice.NevtSpectra(DM_flux_STD,1,years,sparam)
                Emu             = map(lambda EE : ice.MuonEnergy(EE*sparam.GeV,1,sparam)/sparam.GeV, Enu)
            else :
                print "Wrong neutrino type."
                quit()
                
            # converting to "italian" scale
            #evt_spectra_logmuE_STD = []
            evt_spectra_muE_STD = []
            for i,E in enumerate(Enu):
                #evt_spectra_logmuE_STD.append(np.log(10.0)*E*evt_spectra_STD[i])
                evt_spectra_muE_STD.append((E/Emu[i])*evt_spectra_STD[i])
                
            # plot
            #plt.plot(Emu,evt_spectra_logmuE_STD,label = chname[chid], linestyle = "solid", color = colors[chid])
            plt.plot(Emu,evt_spectra_muE_STD,label = chname[chid]+" "+param.name, linestyle = "dashed", color = colors[chid])
        
        plt.title(r"$\chi_m = "+str(DMm/sparam.GeV)+" \mathrm{GeV}$")
        #plt.ylabel(r"$\frac{\mathrm{dN_{\mu}}}{d\mathrm{log_{10}}E_\mu}$")
        plt.ylabel(r"$\frac{\mathrm{dN_{\mu}}}{d E_\mu} [\mathrm{GeV}^{-1}]$")
        plt.xlabel(r"$E_\mu [\mathrm{GeV}]$")
        plt.semilogx()
        plt.xlim(3.0,DMm/sparam.GeV)
        
        #DM_flux_STD, evt_spectra_STD = 0,0
        
    
    fig.subplots_adjust(left=0.1, right=0.85,wspace = 0.35, top = 0.85, bottom = 0.15)
    mpl.rcParams['legend.fontsize'] = "small"
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,fancybox = True)
    
    plt.suptitle("Muon event spectra for an Icecube-like detector")
    
    path = "../plots/"
    #filename = "PlotDM_Icecube_muons_italian_scale.png"
    filename = "PlotDM_Icecube_muons.png"
    plt.savefig(path+filename)
    
    
def PlotDM_MTonDetector(param,sparam = PC.PhysicsConstants()):
    """ Plot muon event spectra from DM annihilation at the Sun on selected
    channels. Can compare to non-standard neutrino oscillation scenarios.
    Uses contained event in a megaton detector.
    
    @type   param   :   physicsconstants
    @param  param   :   set of physical parameters used to make the plot.
    @type   sparam  :   physicsconstants
    @param  sparam  :   standard parameters
    
    @rtype          :   plot
    @return         :   generates the spectra plot
    """
    Sun        = bd.Sun()
    VC         = bd.Vacuum()
    years       = 1.0
    
    DM_mass     = [50.0*sparam.GeV,100.0*sparam.GeV,1.0*sparam.TeV]  # DARK MATTER MASS
    channels    = {'bb' : 0, 'tautau' : 1}#, 'cc': 2}
    chname      = ['$b\\bar{b}$','$\\tau\\bar{\\tau}$','$c\\bar{c}$']
    colors      = ['orange', 'r', 'k', 'c', 'm', 'y', 'k']
    
    fig = plt.figure(figsize = (3*len(DM_mass)+2,4))
    
    for i,DMm in enumerate(DM_mass):
        fig.add_subplot(1,len(DM_mass),i+1)
        for ch,chid in channels.iteritems():
            # define muon neutrinos flux from DM annhiliation
            Enu                 = np.arange(1.0,DMm/param.GeV,0.1)
            
            ##### STANDARD PARAMETER #####
            
            sparam.Refresh()
            
            # calculate muon event spectra at icecube-like detector
            evt_spectra_STD = map(lambda EE : DM.DMneuDet(1,EE*sparam.GeV,ch,DMm,Sun,sparam,True)*(sparam.GeV/DMm),Enu)
            
            if sparam.neutype == "neutrino":
                Emu             = map(lambda EE : ice.MuonEnergy(EE*sparam.GeV,0,sparam)/sparam.GeV, Enu)
            elif sparam.neutype == "antineutrino":
                Emu             = map(lambda EE : ice.MuonEnergy(EE*sparam.GeV,1,sparam)/sparam.GeV, Enu)
            else :
                print "Wrong neutrino type."
                quit()
                
            # converting to "italian" scale
            #evt_spectra_logmuE_STD = []
            evt_spectra_muE_STD = []
            for i,E in enumerate(Enu):
                #evt_spectra_logmuE_STD.append(np.log(10.0)*E*evt_spectra_STD[i])
                evt_spectra_muE_STD.append((E/Emu[i])*evt_spectra_STD[i])
                
            # plot
            #plt.plot(Emu,evt_spectra_logmuE_STD,label = chname[chid], linestyle = "solid", color = colors[chid])
            plt.plot(Emu,evt_spectra_muE_STD,label = chname[chid]+" "+sparam.name, linestyle = "solid", color = colors[chid])

            
            ##### TEST PARAMETER #####
            
            param.Refresh()
            
            # calculate muon event spectra at icecube-like detector
            evt_spectra_STD = map(lambda EE : DM.DMneuDet(1,EE*sparam.GeV,ch,DMm,Sun,param,True)*(sparam.GeV/DMm),Enu)
            
            if sparam.neutype == "neutrino":
                Emu             = map(lambda EE : ice.MuonEnergy(EE*sparam.GeV,0,sparam)/sparam.GeV, Enu)
            elif sparam.neutype == "antineutrino":
                Emu             = map(lambda EE : ice.MuonEnergy(EE*sparam.GeV,1,sparam)/sparam.GeV, Enu)
            else :
                print "Wrong neutrino type."
                quit()
                
            # converting to "italian" scale
            #evt_spectra_logmuE_STD = []
            evt_spectra_muE_STD = []
            for i,E in enumerate(Enu):
                #evt_spectra_logmuE_STD.append(np.log(10.0)*E*evt_spectra_STD[i])
                evt_spectra_muE_STD.append((E/Emu[i])*evt_spectra_STD[i])
                
            # plot
            #plt.plot(Emu,evt_spectra_logmuE_STD,label = chname[chid], linestyle = "solid", color = colors[chid])
            plt.plot(Emu,evt_spectra_muE_STD,label = chname[chid]+" "+param.name  , linestyle = "dashed", color = colors[chid])
            
            
            ##### VACUUM #####
            
            param.Refresh()
            
            # calculate muon event spectra at icecube-like detector
            evt_spectra_VC = map(lambda EE : DM.DMneuDet(1,EE*sparam.GeV,ch,DMm,VC,param,True)*(sparam.GeV/DMm),Enu)
            
            if sparam.neutype == "neutrino":
                Emu             = map(lambda EE : ice.MuonEnergy(EE*sparam.GeV,0,sparam)/sparam.GeV, Enu)
            elif sparam.neutype == "antineutrino":
                Emu             = map(lambda EE : ice.MuonEnergy(EE*sparam.GeV,1,sparam)/sparam.GeV, Enu)
            else :
                print "Wrong neutrino type."
                quit()
                
            # converting to "italian" scale
            #evt_spectra_logmuE_STD = []
            evt_spectra_muE_VC = []
            for i,E in enumerate(Enu):
                #evt_spectra_logmuE_STD.append(np.log(10.0)*E*evt_spectra_STD[i])
                evt_spectra_muE_VC.append((E/Emu[i])*evt_spectra_VC[i])
                
            # plot
            #plt.plot(Emu,evt_spectra_logmuE_STD,label = chname[chid], linestyle = "solid", color = colors[chid])
            plt.plot(Emu,evt_spectra_muE_VC,label = chname[chid]+" "+param.name+" Vacuum"  , linestyle = "dotted", color = colors[chid])
            
            
        
        plt.title(r"$\chi_m = "+str(DMm/sparam.GeV)+" \mathrm{GeV}$")
        #plt.ylabel(r"$\frac{\mathrm{dN_{\mu}}}{d\mathrm{log_{10}}E_\mu}$")
        plt.ylabel(r"$\frac{\mathrm{dN_{\mu}}}{d E_\mu} [\mathrm{GeV}^{-1}]$")
        plt.xlabel(r"$E_\mu [\mathrm{GeV}]$")
        plt.semilogx()
        plt.xlim(3.0,DMm/sparam.GeV)
        
        #DM_flux_STD, evt_spectra_STD = 0,0
        
    
    fig.subplots_adjust(left=0.1, right=0.85,wspace = 0.35, top = 0.85, bottom = 0.15)
    mpl.rcParams['legend.fontsize'] = "small"
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,fancybox = True)
    
    plt.suptitle("Muon event spectra for detector with a Mton-year exposure")
    
    path = "../plots/"
    #filename = "PlotDM_Icecube_muons_italian_scale.png"
    filename = "PlotDM_megaton-year_det_muons_"+sparam.neutype+".png"
    plt.savefig(path+filename)
    
    
#===============================================================================
# MUON EVENT SPECTRA AT ICECUBE-LIKE DETECTOR USING MC DATA
#===============================================================================


def PlotDMMC_Icecube(param,sparam = PC.PhysicsConstants()):
    """ Plot muon event spectra from DM annihilation at the Sun on selected
    channels. Can compare to non-standard neutrino oscillation scenarios.
    Uses Icecube simulation on L{ice.NevtSpectra}. Considering a 1 year
    exposure time. Using MC spectra.
    
    @type   param   :   physicsconstants
    @param  param   :   set of physical parameters used to make the plot.
    @type   sparam  :   physicsconstants
    @param  sparam  :   standard parameters
    
    @rtype          :   plot
    @return         :   generates the spectra plot
    """
    years       = 1.0
    
    DM_mass     = [50.0*sparam.GeV,100.0*sparam.GeV,1.0*sparam.TeV]  # DARK MATTER MASS
    DMsig       = 1.0e-36*param.cm**2
    channels    = {'bb' : 0, 'tautau' : 1}#, 'cc': 2}
    chname      = ['$b\\bar{b}$','$\\tau\\bar{\\tau}$','$c\\bar{c}$']
    colors      = ['orange', 'r', 'k', 'c', 'm', 'y', 'k']
    
    fig = plt.figure(figsize = (3*len(DM_mass)+2,4))
    
    for i,DMm in enumerate(DM_mass):
        fig.add_subplot(1,len(DM_mass),i+1)
        for ch,chid in channels.iteritems():
            
            ## STANRDARD PARAMETER SET
            sparam.Refresh()
            
            # define muon neutrinos flux from DM annhiliation
            DM_flux_STD         = DM.DMFluxAtDetector(ch,DMm,DMsig,sparam)
            Enu                 = np.arange(DM_flux_STD.minenergy/param.GeV,DM_flux_STD.maxenergy/param.GeV,0.1)
            
            # calculate muon event spectra at icecube-like detector
            if sparam.neutype == "neutrino":
                evt_spectra_STD = ice.NevtSpectra(DM_flux_STD,0,years,sparam)
                Emu             = map(lambda EE : ice.MuonEnergy(EE*sparam.GeV,0,sparam)/sparam.GeV, Enu)
            elif sparam.neutype == "antineutrino":
                evt_spectra_STD = ice.NevtSpectra(DM_flux_STD,1,years,sparam)
                Emu             = map(lambda EE : ice.MuonEnergy(EE*sparam.GeV,1,sparam)/sparam.GeV, Enu)
            else :
                print "Wrong neutrino type."
                quit()
                
            # converting to "italian" scale
            #evt_spectra_logmuE_STD = []
            evt_spectra_muE_STD = []
            for i,E in enumerate(Enu):
                #evt_spectra_logmuE_STD.append(np.log(10.0)*E*evt_spectra_STD[i])
                evt_spectra_muE_STD.append(float((E/Emu[i])*evt_spectra_STD[i]))
            
            #print evt_spectra_muE_STD
                                                    
            # plot
            #plt.plot(Emu,evt_spectra_logmuE_STD,label = chname[chid], linestyle = "solid", color = colors[chid])
            plt.plot(Emu,evt_spectra_muE_STD,label = chname[chid]+" "+sparam.name, linestyle = "solid", color = colors[chid])
            
            ## OTHER PARAMETER SET
            param.Refresh()
            
            # define muon neutrinos flux from DM annhiliation
            #DM_flux_STD         = DM.DMFluxAtDetector(ch,DMm,param)
            DM_flux_STD_MC      = DM.DMMCFluxAtDetector(ch,DMm,DMsig,param)[0]
            Enu                 = np.arange(DM_flux_STD_MC.minenergy/sparam.GeV,DM_flux_STD_MC.maxenergy/sparam.GeV,0.1)
            
            # calculate muon event spectra at icecube-like detector
            if sparam.neutype == "neutrino":
                evt_spectra_STD = ice.NevtSpectra(DM_flux_STD_MC,0,years,sparam)
                Emu             = map(lambda EE : ice.MuonEnergy(EE*sparam.GeV,0,sparam)/sparam.GeV, Enu)
            elif sparam.neutype == "antineutrino":
                evt_spectra_STD = ice.NevtSpectra(DM_flux_STD_MC,1,years,sparam)
                Emu             = map(lambda EE : ice.MuonEnergy(EE*sparam.GeV,1,sparam)/sparam.GeV, Enu)
            else :
                print "Wrong neutrino type."
                quit()
                
            # converting to "italian" scale
            #evt_spectra_logmuE_STD = []
            evt_spectra_muE_STD = []
            for i,E in enumerate(Enu):
                #evt_spectra_logmuE_STD.append(np.log(10.0)*E*evt_spectra_STD[i])
                evt_spectra_muE_STD.append((E/Emu[i])*evt_spectra_STD[i])
                
            # plot
            #plt.plot(Emu,evt_spectra_logmuE_STD,label = chname[chid], linestyle = "solid", color = colors[chid])
            plt.plot(Emu,evt_spectra_muE_STD,label = chname[chid]+" "+param.name, linestyle = "dashed", color = colors[chid])
        
        plt.title(r"$\chi_m = "+str(DMm/sparam.GeV)+" \mathrm{GeV}$")
        #plt.ylabel(r"$\frac{\mathrm{dN_{\mu}}}{d\mathrm{log_{10}}E_\mu}$")
        plt.ylabel(r"$\mathrm{Events}$")
        plt.xlabel(r"$E_\mu [\mathrm{GeV}]$")
        plt.semilogx()
        plt.xlim(3.0,DMm/sparam.GeV)
        
        #DM_flux_STD, evt_spectra_STD = 0,0
        
    
    fig.subplots_adjust(left=0.1, right=0.85,wspace = 0.35, top = 0.85, bottom = 0.15)
    mpl.rcParams['legend.fontsize'] = "small"
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,fancybox = True)
    
    plt.suptitle("Muon event spectra for an Icecube-like detector")
    
    path = "../plots/"
    #filename = "PlotDM_Icecube_muons_italian_scale.png"
    filename = "PlotDM_Icecube_muons_MC.png"
    plt.savefig(path+filename)
    
    
def PlotDMMC_MTonDetector(param,sparam = PC.PhysicsConstants()):
    """ Plot muon event spectra from DM annihilation at the Sun on selected
    channels. Can compare to non-standard neutrino oscillation scenarios.
    Uses contained event in a megaton detector. Using MC spectra.
    
    @type   param   :   physicsconstants
    @param  param   :   set of physical parameters used to make the plot.
    @type   sparam  :   physicsconstants
    @param  sparam  :   standard parameters
    
    @rtype          :   plot
    @return         :   generates the spectra plot
    """
    Sun        = bd.Sun()
    VC         = bd.Vacuum()
    years      = 1.0
    DMsig      = 1.0e-40*param.cm**2
    
    
    #DM_mass     = [50.0*sparam.GeV,100.0*sparam.GeV,500.0*sparam.GeV,1.0*sparam.TeV]  # DARK MATTER MASS
    DM_mass     = [500.0*sparam.GeV,1.0*sparam.TeV]  # DARK MATTER MASS
    DM_mass     = [1.0*sparam.TeV]  # DARK MATTER MASS
    channels    = {'bb' : 0, 'tautau' : 1}#, 'cc': 2}
    chname      = ['$b\\bar{b}$','$\\tau\\bar{\\tau}$','$c\\bar{c}$']
    colors      = ['orange', 'r', 'k', 'c', 'm', 'y', 'k']
    
    fig = plt.figure(figsize = (3*len(DM_mass)+2,4))
    
    E_nu_list_1T  = gt.LogSpaceEnergies(1.0,1000.0,binnum = 30)
    E_nu_hpl_1T   = gt.MidPoint(gt.LogSpaceEnergies(1.0,1000.0,binnum = 30))     
    
    for i,DMm in enumerate(DM_mass):
        fig.add_subplot(1,len(DM_mass),i+1)
        for ch,chid in channels.iteritems():
            # define muon neutrinos flux from DM annhiliation
            Enu                 = np.arange(1.0,DMm/param.GeV,0.1)
            E_nu_list  = filter(lambda x : x<DMm/param.GeV, E_nu_list_1T)
            E_nu_hpl   = filter(lambda x : x<=DMm/param.GeV, E_nu_hpl_1T)            
            
            ##### STANDARD PARAMETER #####
            
            #######sparam.Refresh()
            
            ######## calculate muon event spectra at icecube-like detector
            ######## converting to event number. "Eventeando".
            #######
            #######evt_spectra_STD = map(lambda EE : DM.DMneuDet(1,EE*sparam.GeV,ch,DMm,DMsig,Sun,sparam,True)*(sparam.GeV/DMm),Enu)
            #######
            #######Enu_aux = np.arange(1.0,DMm/param.GeV+0.1,0.1)
            #######evt_spectra_STD.append(0.0)
            #######
            #######inter_N_dE = interpolate.interp1d(Enu_aux,evt_spectra_STD)
            #######E_nu_list  = filter(lambda x : x<DMm/param.GeV, E_nu_list_1T)
            #######E_nu_hpl   = filter(lambda x : x<=DMm/param.GeV, E_nu_hpl_1T)
            #######
            #######if np.abs(E_nu_list[-1]-DMm/param.GeV) > 0.5:
            #######    E_nu_list.append(DMm/param.GeV)
            #######
            #######E_nu_list[-1] = DMm/param.GeV
            #######
            #######evt_number_STD = []
            #######
            #######for i in range(len(E_nu_list)-1):
            #######    evt_number_STD.append(integrate.quad(inter_N_dE,E_nu_list[i],E_nu_list[i+1])[0])    
            #######
            #######if sparam.neutype == "neutrino":
            #######    #Emu             = map(lambda EE : ice.MuonEnergy(EE*spara9m.GeV,0,sparam)/sparam.GeV, Enu)
            #######    Emu             = map(lambda EE : ice.MuonEnergy(EE*sparam.GeV,0,sparam)/sparam.GeV, E_nu_hpl)                
            #######elif sparam.neutype == "antineutrino":
            #######    #Emu             = map(lambda EE : ice.MuonEnergy(EE*sparam.GeV,1,sparam)/sparam.GeV, Enu)
            #######    Emu             = map(lambda EE : ice.MuonEnergy(EE*sparam.GeV,0,sparam)/sparam.GeV, E_nu_hpl)                
            #######else :
            #######    print "Wrong neutrino type."
            #######    quit()
            #######
            #######plt.plot(E_nu_hpl,evt_number_STD,label = chname[chid]+" "+sparam.name+" ANA", linestyle = "solid", color = colors[chid])
            #######
            ############ TEST PARAMETER #####
            
            ############ BEGIN 2+3 #####
            
            param.Refresh()
            
            param.name = "2+3"
            
            # calculate muon flux at icecube-like detector
            datapath = "../data/myMC/"+param.name+"/"+param.neutype+"/"
            inters = DM.DMFNeuFluxMCDetv2(ch,DMm,DMsig,param,onlyosc = False,datapath = datapath,crosscheck = False)
            
            param.neutype = "neutrino"
            inter_flux_neu  = inters[0]
            evt_spectra_STD_neu  = map(lambda EE : DM.DMneuDetMC(1,EE*sparam.GeV,ch,DMm,DMsig,Sun,param,inter_flux_neu),E_nu_hpl)
            
            param.neutype = "antineutrino"
            inter_flux_aneu = inters[1]
            evt_spectra_STD_aneu = map(lambda EE : DM.DMneuDetMC(1,EE*sparam.GeV,ch,DMm,DMsig,Sun,param,inter_flux_aneu),E_nu_hpl)
            
            Emu_neu  = map(lambda EE : ice.MuonEnergy(EE*sparam.GeV,0,sparam)/sparam.GeV, E_nu_hpl)
            Emu_aneu = map(lambda EE : ice.MuonEnergy(EE*sparam.GeV,1,sparam)/sparam.GeV, E_nu_hpl)
                
            # converting to "italian" scale
            #evt_spectra_logmuE_STD = []
            evt_spectra_muE_STD = []
            evt_spectra_amuE_STD = []
            for i,E in enumerate(E_nu_hpl):
                evt_spectra_muE_STD.append(float(evt_spectra_STD_neu[i]))
                evt_spectra_amuE_STD.append(float(evt_spectra_STD_aneu[i]))
                
            # calculating total muon events
            #inter_mu  = interpolate.interp1d(Emu_neu,evt_spectra_muE_STD)
            #inter_amu = interpolate.interp1d(Emu_aneu,evt_spectra_amuE_STD)
            
            inter_mu  = interpolate.interp1d(E_nu_hpl,evt_spectra_muE_STD)
            inter_amu = interpolate.interp1d(E_nu_hpl,evt_spectra_amuE_STD)
            
            evt_spectra_total_mu =[]
            #for E in Emu_neu :
            for E in E_nu_hpl:
                if E >= Emu_aneu[0]:
                    evt_spectra_total_mu.append(inter_mu(E)+inter_amu(E))
                else :
                    evt_spectra_total_mu.append(inter_mu(E))
                                    
            # plot
            plt.plot(E_nu_hpl,evt_spectra_total_mu,label = chname[chid]+" "+param.name+""  , linestyle = "dashed", color = colors[chid])
            
            param.neutype = "neutrino"            
            
            ############ END 2+3 #####
            
            sparam.Refresh()
            sdatapath = "../data/myMC/"+sparam.name+"/"+sparam.neutype+"/"
            # calculate muon flux at icecube-like detector
            inters = DM.DMFNeuFluxMCDetv2(ch,DMm,DMsig,sparam,onlyosc = False,datapath = sdatapath)
            
            sparam.neutype = "neutrino"
            inter_flux_neu  = inters[0]
            evt_spectra_STD_neu  = map(lambda EE : DM.DMneuDetMC(1,EE*sparam.GeV,ch,DMm,DMsig,Sun,sparam,inter_flux_neu),E_nu_hpl)
            
            sparam.neutype = "antineutrino"
            inter_flux_aneu = inters[1]
            evt_spectra_STD_aneu = map(lambda EE : DM.DMneuDetMC(1,EE*sparam.GeV,ch,DMm,DMsig,Sun,sparam,inter_flux_aneu),E_nu_hpl)
            
            Emu_neu  = map(lambda EE : ice.MuonEnergy(EE*sparam.GeV,0,sparam)/sparam.GeV, E_nu_hpl)
            Emu_aneu = map(lambda EE : ice.MuonEnergy(EE*sparam.GeV,1,sparam)/sparam.GeV, E_nu_hpl)
                
            # converting to "italian" scale
            #evt_spectra_logmuE_STD = []
            evt_spectra_muE_STD = []
            evt_spectra_amuE_STD = []
            for i,E in enumerate(E_nu_hpl):
                evt_spectra_muE_STD.append(float(evt_spectra_STD_neu[i]))
                evt_spectra_amuE_STD.append(float(evt_spectra_STD_aneu[i]))
                
            # calculating total muon events
            #inter_mu  = interpolate.interp1d(Emu_neu,evt_spectra_muE_STD)
            #inter_amu = interpolate.interp1d(Emu_aneu,evt_spectra_amuE_STD)
            
            inter_mu  = interpolate.interp1d(E_nu_hpl,evt_spectra_muE_STD)
            inter_amu = interpolate.interp1d(E_nu_hpl,evt_spectra_amuE_STD)
            
            evt_spectra_total_mu =[]
            #for E in Emu_neu :
            for E in E_nu_hpl:
                if E >= Emu_aneu[0]:
                    evt_spectra_total_mu.append(inter_mu(E)+inter_amu(E))
                else :
                    evt_spectra_total_mu.append(inter_mu(E))
                                    
            # plot
            plt.plot(E_nu_hpl,evt_spectra_total_mu,label = chname[chid]+" "+param.name+""  , linestyle = "solid", color = colors[chid])
            
            sparam.neutype = "neutrino"
            
            ###### ONLY OSC MC #####
            #
            ## calculate muon flux at icecube-like detector
            #inters = DM.DMFNeuFluxMCDetv2(ch,DMm,DMsig,sparam,onlyosc = True)
            #
            #sparam.neutype = "neutrino"
            #inter_flux_neu  = inters[0]
            #evt_spectra_STD_neu  = map(lambda EE : DM.DMneuDetMC(1,EE*sparam.GeV,ch,DMm,DMsig,Sun,sparam,inter_flux_neu),E_nu_hpl)
            #
            #sparam.neutype = "antineutrino"
            #inter_flux_aneu = inters[1]
            #evt_spectra_STD_aneu = map(lambda EE : DM.DMneuDetMC(1,EE*sparam.GeV,ch,DMm,DMsig,Sun,sparam,inter_flux_aneu),E_nu_hpl)
            #
            #Emu_neu  = map(lambda EE : ice.MuonEnergy(EE*sparam.GeV,0,sparam)/sparam.GeV, E_nu_hpl)
            #Emu_aneu = map(lambda EE : ice.MuonEnergy(EE*sparam.GeV,1,sparam)/sparam.GeV, E_nu_hpl)
            #    
            ## converting to "italian" scale
            ##evt_spectra_logmuE_STD = []
            #evt_spectra_muE_STD = []
            #evt_spectra_amuE_STD = []
            #for i,E in enumerate(E_nu_hpl):
            #    evt_spectra_muE_STD.append(float(evt_spectra_STD_neu[i]))
            #    evt_spectra_amuE_STD.append(float(evt_spectra_STD_aneu[i]))
            #
            #inter_mu  = interpolate.interp1d(E_nu_hpl,evt_spectra_muE_STD)
            #inter_amu = interpolate.interp1d(E_nu_hpl,evt_spectra_amuE_STD)
            #
            #evt_spectra_total_mu =[]
            ##for E in Emu_neu :
            #for E in E_nu_hpl:
            #    if E >= Emu_aneu[0]:
            #        evt_spectra_total_mu.append(inter_mu(E)+inter_amu(E))
            #    else :
            #        evt_spectra_total_mu.append(inter_mu(E))
            #                        
            ## plot
            #plt.plot(E_nu_hpl,evt_spectra_total_mu,label = chname[chid]+" "+param.name+" MC OO"  , linestyle = "solid", color = colors[chid])
            #
            #sparam.neutype = "neutrino"            
            #
            
        plt.title(r"$\chi_m = "+str(DMm/sparam.GeV)+" \mathrm{GeV}$")
        plt.ylabel(r"$Events$")
        plt.xlabel(r"$E_\nu [\mathrm{GeV}]$")
        plt.semilogx()
        plt.xlim(3.0,DMm/sparam.GeV)
        
        #DM_flux_STD, evt_spectra_STD = 0,0
        
    
    fig.subplots_adjust(left=0.1, right=0.85,wspace = 0.35, top = 0.85, bottom = 0.15)
    mpl.rcParams['legend.fontsize'] = "small"
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,fancybox = True)
    
    plt.suptitle("Muon event spectra for detector with a Mton-year exposure")
    
    path = "../plots/"
    #filename = "PlotDM_Icecube_muons_italian_scale.png"
    filename = "PlotDM_megaton-year_det_muons_events_"+sparam.neutype+"_MC.png"
    plt.savefig(path+filename)
    
def PlotFluxAtDetector(param,sparam = PC.PhysicsConstants()):
    """ Plot muon event spectra from DM annihilation at the Sun on selected
    channels. Can compare to non-standard neutrino oscillation scenarios.
    Uses contained event in a megaton detector. Using MC spectra.
    
    @type   param   :   physicsconstants
    @param  param   :   set of physical parameters used to make the plot.
    @type   sparam  :   physicsconstants
    @param  sparam  :   standard parameters
    
    @rtype          :   plot
    @return         :   generates the spectra plot
    """
    Sun        = bd.Sun()
    DMsig      = 1.0e-40*param.cm**2
    
    DM_mass     = [50.0*sparam.GeV,100.0*sparam.GeV,1.0*sparam.TeV]  # DARK MATTER MASS
    channels    = {'bb' : 0, 'tautau' : 1}#, 'cc': 2}
    chname      = ['$b\\bar{b}$','$\\tau\\bar{\\tau}$','$c\\bar{c}$']
    colors      = ['orange', 'r', 'k', 'c', 'm', 'y', 'k']
    
    fig = plt.figure(figsize = (3*len(DM_mass)+2,4))
    
    E_nu_list_1T  = gt.LogSpaceEnergies(1.0,1000.0,binnum = 30)
    E_nu_hpl_1T   = gt.MidPoint(gt.LogSpaceEnergies(1.0,1000.0,binnum = 30))     
    
    
    for i,DMm in enumerate(DM_mass):
        fig.add_subplot(1,len(DM_mass),i+1)
        for ch,chid in channels.iteritems():
            Enu                 = np.arange(1.0,DMm/param.GeV,0.1)
            
            ##### STANDARD PARAMETER #####
            
            sparam.Refresh()
            
            # calculate muon event spectra at icecube-like detector
            
            flux_spectra_STD = map(lambda EE : DM.DMFluxneuDet(1,EE*sparam.GeV,ch,DMm,DMsig,Sun,sparam,True)*(sparam.GeV/DMm),Enu)
            
            Enu_aux = np.arange(1.0,DMm/param.GeV+0.1,0.1)
            flux_spectra_STD.append(0.0)
            
            inter_dPhi_dE = interpolate.interp1d(Enu_aux,flux_spectra_STD)
            E_nu_list  = filter(lambda x : x<DMm/param.GeV, E_nu_list_1T)
            E_nu_hpl   = filter(lambda x : x<=DMm/param.GeV, E_nu_hpl_1T)
            
            if np.abs(E_nu_list[-1]-DMm/param.GeV) > 0.5:
                E_nu_list.append(DMm/param.GeV)
            
            E_nu_list[-1] = DMm/param.GeV
            
            flux_bin_STD = []
            
            for i in range(len(E_nu_list)-1):
                flux_bin_STD.append(integrate.quad(inter_dPhi_dE,E_nu_list[i],E_nu_list[i+1])[0])

            plt.plot(E_nu_hpl,flux_bin_STD,label = chname[chid]+" "+sparam.name, linestyle = "solid", color = colors[chid])
            #plt.plot(Enu,flux_spectra_STD,label = chname[chid]+" "+sparam.name, linestyle = "solid", color = colors[chid])
            
            ##### TEST PARAMETER #####
            
            param.Refresh()
            
            inters = DM.DMFNeuFluxMCDetv2(ch,DMm,DMsig,param)
            
            inter_flux_neu  = inters[0]
            inter_flux_aneu = inters[1]
                
            flux_spectra = [inter_flux_neu(E) for E in E_nu_hpl]
            #flux_spectra = [inter_flux_neu(E)+inter_flux_aneu(E) for E in E_nu_hpl]
                
            # plot
            plt.plot(E_nu_hpl,flux_spectra,label = chname[chid]+" "+param.name  , linestyle = "dashed", color = colors[chid])
            
            param.neutype = "neutrino"
            
        plt.title(r"$\chi_m = "+str(DMm/sparam.GeV)+" \mathrm{GeV}$")
        plt.ylabel(r"$\frac{d\phi}{dE_\nu}$")
        plt.xlabel(r"$E_\nu [\mathrm{GeV}]$")
        plt.semilogx()
        plt.xlim(3.0,DMm/sparam.GeV)
    
    fig.subplots_adjust(left=0.1, right=0.85,wspace = 0.35, top = 0.85, bottom = 0.15)
    mpl.rcParams['legend.fontsize'] = "small"
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,fancybox = True)
    
    plt.suptitle(r"$\nu_\mu$ flux at detector")
    
    path = "../plots/"
    filename = "PlotDM_megaton-year_det_flux_MC.png"
    plt.savefig(path+filename)
    
def PlotNewIcecubeLimit(param):
    """Plots new Icecube limit considering 2+3 case on the flux
    # param         : physics parameter
    """
    DMm     = [250.0*param.GeV,500.0*param.GeV,1000.0*param.GeV]
    DMsig_soft   = [0.0,2.5e-41*param.cm**2,1.8e-41*param.cm**2]
    DMsig_hard   = [3.70e-43*param.cm**2,2.9e-43*param.cm**2,7.2e-43*param.cm**2]
    
    #DMm = [500.0*param.GeV]
    #DMsig_hard = [2.9e-43*param.cm**2]
    
    savedatapath = "../data/DMConvolutedFlux/"
    
    channels = ['bb','WW']
    
    glim = []
    flim = []
    
    for ch in channels : 
        if ch == 'bb':
            DMsig = DMsig_soft
        elif ch == 'WW':
            DMsig = DMsig_hard
        else :
            print "Missing cross sections for this channel."
            quit()
            
        dat_pair = [[DMm[i],DMsig[i]] for i in range(len(DMm))]        
        
        E_nu_list_1T  = gt.LogSpaceEnergies(1.0,1000.0,binnum = 30)
        E_nu_hpl_1T   = gt.MidPoint(gt.LogSpaceEnergies(1.0,1000.0,binnum = 30))    
        
        gamma_lim = []
        flux_lim  = []
        
        for dm, sig in dat_pair:
            
            E_nu_list  = filter(lambda x : x<dm/param.GeV, E_nu_list_1T) # [GeV]
            E_nu_hpl   = filter(lambda x : x<=dm/param.GeV, E_nu_hpl_1T) # [GeV]
            
            #print dm,sig
            DM_annihilation_rate_Sun = DM.DMSunAnnihilationRate(dm,sig,param)     # [eV]
            
            print np.sum(DM_annihilation_rate_Sun)*param.sec
            normalization = np.sum((DM_annihilation_rate_Sun/(4.0*np.pi*param.AU**2)))  # [eV^3]
            
            try :
                print "reading data"
                filename_flux  = "DM_flux_"+str(dm/param.GeV)+"_GeV_xsection_"+str(sig/param.cm**2)+"_cm2_osc_"+param.name+"_channel_"+ch+".dat"
                filename_gamma = "DM_gamma_"+str(dm/param.GeV)+"_GeV_xsection_"+str(sig/param.cm**2)+"_cm2_osc_"+param.name+"_channel_"+ch+".dat"
                
                file_flux  = open(savedatapath+filename_flux, 'r')
                file_gamma = open(savedatapath+filename_gamma, 'r')
                
                fluxint  = []
                gammaint = []
                
                gt.hreadfilev4(file_flux,fluxint,param)
                gt.hreadfilev4(file_gamma,gammaint,param)
                
                fluxint = fluxint[0]
                gammaint = gammaint[0]
                
                file_flux.close()
                file_gamma.close()
            except :
                print "generating data"
                # STD
                param.neutype = "neutrino"
                param.name = "STD"            
                ## 2+3 - neutrino
                #param.neutype = "neutrino"
                #param.name = "2+3"
                ## 2+3 - antineutrino
                #param.neutype = "antineutrino"
                #param.name = "2+3"
                
                datapath = "../data/myMC/"+param.name+"/"+param.neutype+"/"
                
                inter_neu,inter_aneu = DM.DMFNeuFluxMCDetv2(ch,dm,sig,param,onlyosc = False, datapath = datapath)
                
                fluxint  = [[E,float(inter_neu(E))*xs.signuNCC(E,0)*param.cm**2*(0.918*param.gr*param.cm**-3)/(939.27*param.MeV)*ice.MuonRange(ice.MuonEnergy(E*param.GeV,0,param),param)*param.meter] for E in E_nu_hpl]
                gammaint = [[E,float(inter_neu(E))*xs.signuNCC(E,0)*param.cm**2*(0.918*param.gr*param.cm**-3)/(939.27*param.MeV)] for E in E_nu_hpl]
                  
                #begin saving flux
                filename_flux  = "DM_flux_"+str(dm/param.GeV)+"_GeV_xsection_"+str(sig/param.cm**2)+"_cm2_osc_"+param.name+"_channel_"+ch+".dat"
                filename_gamma = "DM_gamma_"+str(dm/param.GeV)+"_GeV_xsection_"+str(sig/param.cm**2)+"_cm2_osc_"+param.name+"_channel_"+ch+".dat"
                
                file_flux  = open(savedatapath+filename_flux, 'w')
                file_gamma = open(savedatapath+filename_gamma, 'w')
                
                gt.hwritefile(file_flux ,fluxint ,param)
                gt.hwritefile(file_gamma,gammaint,param)            
                #end saving flux
                file_flux.close()
                file_gamma.close()
            
            pflux  = [f[1] for f in fluxint]
            pgamma = [f[1] for f in gammaint]
            
            int_flux  = sum(pflux)
            int_gamma = sum(pgamma)
            
            flux_lim.append(int_flux*param.km**2*param.year)
            gamma_lim.append(int_gamma*param.km**3*param.year)
            
        flim.append(flux_lim)
        glim.append(gamma_lim)
        
    print flim
    print glim
    
    # Making new limits plot
    
    my_nu_flux_lim_soft = flim[0]
    my_nu_flux_lim_hard = flim[1]    
    
    # from Ref.: arXiv 0902.2460
    DM_mass     = [250.0,500.0,1000.0]
    nu_flux_lim_hard = [8.8e2, 4.2e2, 3.6e2] # km^-2y^-1
    nu_flux_lim_soft = [0.0  , 3.5e3, 1.3e3] # km^-2y^-1
    
    # my limits
    
    fig = plt.figure()
    
    plt.plot(DM_mass,nu_flux_lim_hard,'s',linestyle = "solid" ,label = "Icecube (hard)", color = "red")
    plt.plot(DM_mass,nu_flux_lim_soft,'s',linestyle = "dashed",label = "Icecube (soft)", color = "red")
    
    plt.plot(DM_mass,my_nu_flux_lim_hard,'o',linestyle = "solid" ,label = "My Limit (hard)", color = "blue")
    plt.plot(DM_mass,my_nu_flux_lim_soft,'o',linestyle = "dashed",label = "My Limit (soft)", color = "blue")    

    plt.loglog()
    
    plt.xlim(10.0,1.0e4)
    plt.ylim(50.0,1.0e5)
    
    plt.ylabel('$\mathrm{Muon\ flux\ from\ the\ Sun\ [km^{-2}y^{-1}]}$')
    plt.xlabel('$\mathrm{Neutralino\ mass\ [GeV]}$')
    
    plt.legend(loc = 'upper left')
    
    path = "../plots/"
    filename = "PlotNewICFluxLimit.png"
    
    plt.savefig(path+filename)
    
def PlotNewIcecubeXsectionLimit(param, sparam = PC.PhysicsConstants(),param2 = None, use_MC_data = True):
    """Plots new Icecube limit considering 2+3 case on the flux
    # param         : physics parameter
    # use_MC_data   : Toogle the use of MC data. Else will use analitic aproximation.
    """
    #DMm     = [250.0*param.GeV,500.0*param.GeV,1000.0*param.GeV]
    DMm     = [250.0*param.GeV,500.0*param.GeV,1000.0*param.GeV,3000.0*param.GeV,5000.0*param.GeV]
    DMsig_soft   = [0.0,2.5e-41*param.cm**2,1.8e-41*param.cm**2,5.3e-41*param.cm**2,1.1e-40*param.cm**2]
    DMsig_hard   = [3.70e-43*param.cm**2,2.9e-43*param.cm**2,7.2e-43*param.cm**2,7.4e-42*param.cm**2,2.6e-41*param.cm**2]

    savedatapath = "../data/DMConvolutedFlux/"
    
    channels = ['bb','WW']
    
    glim = []
    flim = []
    
    sglim = []
    sflim = []    
    
    p2flim = []
    p2glim = []
    
    for ch in channels : 
        if ch == 'bb':
            DMsig = DMsig_soft
        elif ch == 'WW':
            DMsig = DMsig_hard
        else :
            print "Missing cross sections for this channel."
            quit()
            
        dat_pair = [[DMm[i],DMsig[i]] for i in range(len(DMm))]        
        
        binnum = 100
        use_old_data = True
        
        E_nu_list_1T  = gt.LogSpaceEnergies(1.0,100000.0,binnum = binnum)
        E_nu_hpl_1T   = gt.MidPoint(gt.LogSpaceEnergies(1.0,100000.0,binnum = binnum))    
        
        ########################## FOR PARAM ###################################
        
        param.Refresh()
        
        gamma_lim = []
        flux_lim  = []
        
        for dm, sig in dat_pair:
            
            E_nu_list  = filter(lambda x : x<dm/param.GeV, E_nu_list_1T) # [GeV]
            E_nu_hpl   = filter(lambda x : x<=dm/param.GeV, E_nu_hpl_1T) # [GeV]
            
            #print dm,sig
            DM_annihilation_rate_Sun = DM.DMSunAnnihilationRate(dm,sig,param)     # [eV]
            
            #normalization = np.sum((DM_annihilation_rate_Sun/(4.0*np.pi*param.AU**2)))  # [eV^3]
            
            
            if use_MC_data : 
                try :
                    print "reading data"
                    filename_flux  = "DM_flux_"+str(dm/param.GeV)+"_GeV_xsection_"+str(sig/param.cm**2)+"_cm2_osc_"+param.name+"_channel_"+ch+".dat"
                    filename_gamma = "DM_gamma_"+str(dm/param.GeV)+"_GeV_xsection_"+str(sig/param.cm**2)+"_cm2_osc_"+param.name+"_channel_"+ch+".dat"
                    
                    file_flux  = open(savedatapath+filename_flux, 'r')
                    file_gamma = open(savedatapath+filename_gamma, 'r')
                    
                    fluxint  = []
                    gammaint = []
                    
                    gt.hreadfilev4(file_flux,fluxint,param)
                    gt.hreadfilev4(file_gamma,gammaint,param)
                    
                    fluxint = fluxint[0]
                    gammaint = gammaint[0]
                    
                    file_flux.close()
                    file_gamma.close()
                except :
                    print "generating data"
             
                    param.neutype = "neutrino"
                    datapath = "../data/myMC/"+param.name+"/"+param.neutype+"/"       
                    
                    inter_neu,inter_aneu = DM.DMFNeuFluxMCDetv2(ch,dm,sig,param,use_old_data = use_old_data, datapath = datapath)
                    
                    fluxint_neu  = [[E,float(inter_neu(E))*xs.signuNCC(E,0)*param.cm**2*(0.918*param.gr*param.cm**-3)/(939.27*param.MeV)*ice.MuonRange(ice.MuonEnergy(E*param.GeV,0,param),param)*param.meter
                                     + float(inter_aneu(E))*xs.signuNCC(E,1)*param.cm**2*(0.918*param.gr*param.cm**-3)/(939.27*param.MeV)*ice.MuonRange(ice.MuonEnergy(E*param.GeV,1,param),param)*param.meter
                                    ] for E in E_nu_hpl]
                    gammaint_neu = [[E,float(inter_neu(E))*xs.signuNCC(E,0)*param.cm**2*(0.918*param.gr*param.cm**-3)/(939.27*param.MeV)
                                     + float(inter_neu(E))*xs.signuNCC(E,1)*param.cm**2*(0.918*param.gr*param.cm**-3)/(939.27*param.MeV)
                                    ] for E in E_nu_hpl]                    
                    
                    param.neutype = "antineutrino"
                    datapath = "../data/myMC/"+param.name+"/"+param.neutype+"/"
                    
                    ainter_neu,ainter_aneu = DM.DMFNeuFluxMCDetv2(ch,dm,sig,param,use_old_data = use_old_data, datapath = datapath)                    
                    
                    fluxint_aneu  = [[E,float(ainter_neu(E))*xs.signuNCC(E,0)*param.cm**2*(0.918*param.gr*param.cm**-3)/(939.27*param.MeV)*ice.MuonRange(ice.MuonEnergy(E*param.GeV,0,param),param)*param.meter
                                      + float(ainter_aneu(E))*xs.signuNCC(E,1)*param.cm**2*(0.918*param.gr*param.cm**-3)/(939.27*param.MeV)*ice.MuonRange(ice.MuonEnergy(E*param.GeV,1,param),param)*param.meter
                                    ] for E in E_nu_hpl]
                    gammaint_aneu = [[E,float(ainter_neu(E))*xs.signuNCC(E,0)*param.cm**2*(0.918*param.gr*param.cm**-3)/(939.27*param.MeV)
                                      + float(ainter_aneu(E))*xs.signuNCC(E,1)*param.cm**2*(0.918*param.gr*param.cm**-3)/(939.27*param.MeV)
                                    ] for E in E_nu_hpl]
                    
                    fluxint  = [[E,fluxint_neu[i][1]+fluxint_aneu[i][1]] for i in range(len(E_nu_hpl))]
                    gammaint = [[E,gammaint_neu[i][1]+gammaint_aneu[i][1]] for i in range(len(E_nu_hpl))]
                      
                    #begin saving flux
                    filename_flux  = "DM_flux_"+str(dm/param.GeV)+"_GeV_xsection_"+str(sig/param.cm**2)+"_cm2_osc_"+param.name+"_channel_"+ch+".dat"
                    filename_gamma = "DM_gamma_"+str(dm/param.GeV)+"_GeV_xsection_"+str(sig/param.cm**2)+"_cm2_osc_"+param.name+"_channel_"+ch+".dat"
                    
                    file_flux  = open(savedatapath+filename_flux, 'w')
                    file_gamma = open(savedatapath+filename_gamma, 'w')
                    
                    gt.hwritefile(file_flux ,fluxint ,param)
                    gt.hwritefile(file_gamma,gammaint,param)            
                    #end saving flux
                    file_flux.close()
                    file_gamma.close()
                                        
                pflux  = [f[1] for f in fluxint]
                pgamma = [f[1] for f in gammaint]
                
                int_flux  = sum(pflux)
                int_gamma = sum(pgamma)
                
            else:
                try:
                    print "reading data"
                    filename_flux  = "DM_ANA_flux_"+str(dm/param.GeV)+"_GeV_xsection_"+str(sig/param.cm**2)+"_cm2_osc_"+param.name+"_channel_"+ch+".dat"
                    filename_gamma = "DM_ANA_gamma_"+str(dm/param.GeV)+"_GeV_xsection_"+str(sig/param.cm**2)+"_cm2_osc_"+param.name+"_channel_"+ch+".dat"
                    
                    file_flux  = open(savedatapath+filename_flux, 'r')
                    file_gamma = open(savedatapath+filename_gamma, 'r')
                    
                    fluxint  = []
                    gammaint = []
                    
                    gt.hreadfilev4(file_flux,fluxint,param)
                    gt.hreadfilev4(file_gamma,gammaint,param)
                    
                    fluxint = fluxint[0]
                    gammaint = gammaint[0]
                    
                    file_flux.close()
                    file_gamma.close()
                    
                    int_flux = fluxint[0][0]
                    int_gamma = gammaint[0][0]
                    
                except:
                    print "Data not found. Calculating."
                    # will estimate the fluxes considering propagation without regeneration or NC.
                    datapath = "../data/SunOscProbabilities/test/"
                    mu_inter = DM.DMNeuFluxDetNoInt(ch,dm,sig,param,onlyosc = False, datapath = datapath)
                    # for compatibility
                    int_gamma = 0.0
                    
                    # integrating
                    int_flux = integrate.quad(mu_inter,1.0*param.GeV,dm)[0]
                    
                    gammaint = [int_gamma]
                    fluxint  = [int_flux]
                    
                    #begin saving flux
                    filename_flux  = "DM_ANA_flux_"+str(dm/param.GeV)+"_GeV_xsection_"+str(sig/param.cm**2)+"_cm2_osc_"+param.name+"_channel_"+ch+".dat"
                    filename_gamma = "DM_ANA_gamma_"+str(dm/param.GeV)+"_GeV_xsection_"+str(sig/param.cm**2)+"_cm2_osc_"+param.name+"_channel_"+ch+".dat"
                    
                    file_flux  = open(savedatapath+filename_flux, 'w')
                    file_gamma = open(savedatapath+filename_gamma, 'w')
                    
                    gt.hwritefile(file_flux ,fluxint ,param)
                    gt.hwritefile(file_gamma,gammaint,param)            
                    #end saving flux
                    file_flux.close()
                    file_gamma.close()

            
            flux_lim.append(int_flux*param.km**2*param.year)
            gamma_lim.append(int_gamma*param.km**3*param.year)
            
        flim.append(flux_lim)
        glim.append(gamma_lim)
        
        ########################## FOR SPARAM ###################################
        
        sparam.Refresh()
        
        gamma_lim = []
        flux_lim  = []        
               
        for dm, sig in dat_pair:
            
            E_nu_list  = filter(lambda x : x<dm/param.GeV, E_nu_list_1T) # [GeV]
            E_nu_hpl   = filter(lambda x : x>5.0,filter(lambda x : x<=850.0, E_nu_hpl_1T)) # [GeV]
            
            #print dm,sig
            DM_annihilation_rate_Sun = DM.DMSunAnnihilationRate(dm,sig,sparam)     # [eV]
            
            normalization = np.sum((DM_annihilation_rate_Sun/(4.0*np.pi*param.AU**2)))  # [eV^3]
            
            
            if use_MC_data : 
                try :
                    print "reading data"
                    filename_flux  = "DM_flux_"+str(dm/param.GeV)+"_GeV_xsection_"+str(sig/param.cm**2)+"_cm2_osc_"+sparam.name+"_channel_"+ch+".dat"
                    filename_gamma = "DM_gamma_"+str(dm/param.GeV)+"_GeV_xsection_"+str(sig/param.cm**2)+"_cm2_osc_"+sparam.name+"_channel_"+ch+".dat"
                    
                    file_flux  = open(savedatapath+filename_flux, 'r')
                    file_gamma = open(savedatapath+filename_gamma, 'r')
                    
                    fluxint  = []
                    gammaint = []
                    
                    gt.hreadfilev4(file_flux,fluxint,sparam)
                    gt.hreadfilev4(file_gamma,gammaint,sparam)
                    
                    fluxint = fluxint[0]
                    gammaint = gammaint[0]
                    
                    file_flux.close()
                    file_gamma.close()
                except :
                    print "generating data"
                    # STD
                    sparam.neutype = "neutrino"
                    datapath = "../data/myMC/"+sparam.name+"/"+sparam.neutype+"/"       
                    
                    inter_neu,inter_aneu = DM.DMFNeuFluxMCDetv2(ch,dm,sig,sparam,use_old_data = use_old_data, datapath = datapath)
                    
                    fluxint_neu  = [[E,float(inter_neu(E))*xs.signuNCC(E,0)*sparam.cm**2*(0.918*sparam.gr*sparam.cm**-3)/(939.27*sparam.MeV)*ice.MuonRange(ice.MuonEnergy(E*sparam.GeV,0,sparam),sparam)*sparam.meter
                                     + float(inter_aneu(E))*xs.signuNCC(E,1)*sparam.cm**2*(0.918*sparam.gr*sparam.cm**-3)/(939.27*sparam.MeV)*ice.MuonRange(ice.MuonEnergy(E*sparam.GeV,1,sparam),sparam)*sparam.meter
                                    ] for E in E_nu_hpl]
                    gammaint_neu = [[E,float(inter_neu(E))*xs.signuNCC(E,0)*sparam.cm**2*(0.918*sparam.gr*sparam.cm**-3)/(939.27*sparam.MeV)
                                     + float(inter_aneu(E))*xs.signuNCC(E,1)*sparam.cm**2*(0.918*sparam.gr*sparam.cm**-3)/(939.27*sparam.MeV)
                                    ] for E in E_nu_hpl]                    
                    
                    sparam.neutype = "antineutrino"
                    datapath = "../data/myMC/"+sparam.name+"/"+sparam.neutype+"/"
                    
                    ainter_neu,ainter_aneu = DM.DMFNeuFluxMCDetv2(ch,dm,sig,sparam,use_old_data = use_old_data, datapath = datapath)                    
                    
                    fluxint_aneu  = [[E,float(ainter_neu(E))*xs.signuNCC(E,0)*sparam.cm**2*(0.918*sparam.gr*param.cm**-3)/(939.27*sparam.MeV)*ice.MuonRange(ice.MuonEnergy(E*sparam.GeV,0,sparam),sparam)*sparam.meter
                                      + float(ainter_aneu(E))*xs.signuNCC(E,1)*sparam.cm**2*(0.918*sparam.gr*param.cm**-3)/(939.27*sparam.MeV)*ice.MuonRange(ice.MuonEnergy(E*sparam.GeV,1,sparam),sparam)*sparam.meter
                                    ] for E in E_nu_hpl]
                    gammaint_aneu = [[E,float(ainter_neu(E))*xs.signuNCC(E,0)*sparam.cm**2*(0.918*sparam.gr*param.cm**-3)/(939.27*param.MeV)
                                      + float(ainter_aneu(E))*xs.signuNCC(E,1)*sparam.cm**2*(0.918*sparam.gr*param.cm**-3)/(939.27*param.MeV)
                                    ] for E in E_nu_hpl]
                    
                    fluxint  = [[E,fluxint_neu[i][1]+fluxint_aneu[i][1]] for i in range(len(E_nu_hpl))]
                    gammaint = [[E,gammaint_neu[i][1]+gammaint_aneu[i][1]] for i in range(len(E_nu_hpl))]
                      
                    #begin saving flux
                    filename_flux  = "DM_flux_"+str(dm/sparam.GeV)+"_GeV_xsection_"+str(sig/sparam.cm**2)+"_cm2_osc_"+sparam.name+"_channel_"+ch+".dat"
                    filename_gamma = "DM_gamma_"+str(dm/sparam.GeV)+"_GeV_xsection_"+str(sig/sparam.cm**2)+"_cm2_osc_"+sparam.name+"_channel_"+ch+".dat"
                    
                    file_flux  = open(savedatapath+filename_flux, 'w')
                    file_gamma = open(savedatapath+filename_gamma, 'w')
                    
                    gt.hwritefile(file_flux ,fluxint ,sparam)
                    gt.hwritefile(file_gamma,gammaint,sparam)            
                    #end saving flux
                    file_flux.close()
                    file_gamma.close()
                    
                pflux  = [f[1] for f in fluxint]
                pgamma = [f[1] for f in gammaint]
                
                int_flux  = sum(pflux)
                int_gamma = sum(pgamma)
                
            else:
                try:
                    print "reading data"
                    filename_flux  = "DM_ANA_flux_"+str(dm/param.GeV)+"_GeV_xsection_"+str(sig/param.cm**2)+"_cm2_osc_"+sparam.name+"_channel_"+ch+".dat"
                    filename_gamma = "DM_ANA_gamma_"+str(dm/param.GeV)+"_GeV_xsection_"+str(sig/param.cm**2)+"_cm2_osc_"+sparam.name+"_channel_"+ch+".dat"
                    
                    file_flux  = open(savedatapath+filename_flux, 'r')
                    file_gamma = open(savedatapath+filename_gamma, 'r')
                    
                    fluxint  = []
                    gammaint = []
                    
                    gt.hreadfilev4(file_flux,fluxint,param)
                    gt.hreadfilev4(file_gamma,gammaint,param)
                    
                    fluxint = fluxint[0]
                    gammaint = gammaint[0]
                    
                    file_flux.close()
                    file_gamma.close()
                    
                    int_flux = fluxint[0][0]
                    int_gamma = gammaint[0][0]
                    
                except:
                    print "Data not found. Calculating."
                    # will estimate the fluxes considering propagation without regeneration or NC.
                    datapath = "../data/SunOscProbabilities/"
                    mu_inter = DM.DMNeuFluxDetNoInt(ch,dm,sig,sparam,onlyosc = False, datapath = datapath)
                    
                    # for compatibility
                    int_gamma = 0.0
                    
                    # integrating
                    int_flux = integrate.quad(mu_inter,1.0*param.GeV,dm)[0]
                    
                    gammaint = [int_gamma]
                    fluxint  = [int_flux]
                    
                    #begin saving flux
                    filename_flux  = "DM_ANA_flux_"+str(dm/param.GeV)+"_GeV_xsection_"+str(sig/param.cm**2)+"_cm2_osc_"+sparam.name+"_channel_"+ch+".dat"
                    filename_gamma = "DM_ANA_gamma_"+str(dm/param.GeV)+"_GeV_xsection_"+str(sig/param.cm**2)+"_cm2_osc_"+sparam.name+"_channel_"+ch+".dat"
                    
                    file_flux  = open(savedatapath+filename_flux, 'w')
                    file_gamma = open(savedatapath+filename_gamma, 'w')
                    
                    gt.hwritefile(file_flux ,fluxint ,param)
                    gt.hwritefile(file_gamma,gammaint,param)            
                    #end saving flux
                    file_flux.close()
                    file_gamma.close()                
            
            flux_lim.append(int_flux*param.km**2*param.year)
            gamma_lim.append(int_gamma*param.km**3*param.year)
            
        sflim.append(flux_lim)
        sglim.append(gamma_lim)
        
        ########################## FOR PARAM2 ###################################
        if param2 != None :
            param2.Refresh()
            
            gamma_lim = []
            flux_lim  = []
            
            for dm, sig in dat_pair:
                
                E_nu_list  = filter(lambda x : x<dm/param.GeV, E_nu_list_1T) # [GeV]
                E_nu_hpl   = filter(lambda x : x<=dm/param.GeV, E_nu_hpl_1T) # [GeV]
                
                #print dm,sig
                DM_annihilation_rate_Sun = DM.DMSunAnnihilationRate(dm,sig,param)     # [eV]
                
                #normalization = np.sum((DM_annihilation_rate_Sun/(4.0*np.pi*param.AU**2)))  # [eV^3]
                
                
                if use_MC_data : 
                    try :
                        print "reading data"
                        filename_flux  = "DM_flux_"+str(dm/param.GeV)+"_GeV_xsection_"+str(sig/param.cm**2)+"_cm2_osc_"+param2.name+"_channel_"+ch+".dat"
                        filename_gamma = "DM_gamma_"+str(dm/param.GeV)+"_GeV_xsection_"+str(sig/param.cm**2)+"_cm2_osc_"+param2.name+"_channel_"+ch+".dat"
                        
                        file_flux  = open(savedatapath+filename_flux, 'r')
                        file_gamma = open(savedatapath+filename_gamma, 'r')
                        
                        fluxint  = []
                        gammaint = []
                        
                        gt.hreadfilev4(file_flux,fluxint,param2)
                        gt.hreadfilev4(file_gamma,gammaint,param2)
                        
                        fluxint = fluxint[0]
                        gammaint = gammaint[0]
                        
                        file_flux.close()
                        file_gamma.close()
                    except :
                        print "generating data"
                        # STD
                        #param.neutype = "neutrino"
                        #param.name = "STD"            
                        ## 2+3 - neutrino
                        #param.neutype = "neutrino"
                        #param.name = "2+3"
                        ## 2+3 - antineutrino
                        #param.neutype = "antineutrino"
                        #param.name = "2+3"
                        
                        datapath = "../data/myMC/"+param.name2+"/"+param.neutype+"/"
                        
                        inter_neu,inter_aneu = DM.DMFNeuFluxMCDetv2(ch,dm,sig,param2,use_old_data = True, datapath = datapath)
                        
                        fluxint  = [[E,float(inter_neu(E))*xs.signuNCC(E,0)*param.cm**2*(0.918*param.gr*param.cm**-3)/(939.27*param.MeV)*ice.MuonRange(ice.MuonEnergy(E*param.GeV,0,param2),param2)*param.meter] for E in E_nu_hpl]
                        gammaint = [[E,float(inter_neu(E))*xs.signuNCC(E,0)*param.cm**2*(0.918*param.gr*param.cm**-3)/(939.27*param.MeV)] for E in E_nu_hpl]
                          
                        #begin saving flux
                        filename_flux  = "DM_flux_"+str(dm/param.GeV)+"_GeV_xsection_"+str(sig/param.cm**2)+"_cm2_osc_"+param2.name+"_channel_"+ch+".dat"
                        filename_gamma = "DM_gamma_"+str(dm/param.GeV)+"_GeV_xsection_"+str(sig/param.cm**2)+"_cm2_osc_"+param2.name+"_channel_"+ch+".dat"
                        
                        file_flux  = open(savedatapath+filename_flux, 'w')
                        file_gamma = open(savedatapath+filename_gamma, 'w')
                        
                        gt.hwritefile(file_flux ,fluxint ,param2)
                        gt.hwritefile(file_gamma,gammaint,param2)            
                        #end saving flux
                        file_flux.close()
                        file_gamma.close()
                        
                    pflux  = [f[1] for f in fluxint]
                    pgamma = [f[1] for f in gammaint]
                    
                    int_flux  = sum(pflux)
                    int_gamma = sum(pgamma)
                    
                else:
                    try:
                        print "reading data"
                        filename_flux  = "DM_ANA_flux_"+str(dm/param.GeV)+"_GeV_xsection_"+str(sig/param.cm**2)+"_cm2_osc_"+param2.name+"_channel_"+ch+".dat"
                        filename_gamma = "DM_ANA_gamma_"+str(dm/param.GeV)+"_GeV_xsection_"+str(sig/param.cm**2)+"_cm2_osc_"+param2.name+"_channel_"+ch+".dat"
                        
                        file_flux  = open(savedatapath+filename_flux, 'r')
                        file_gamma = open(savedatapath+filename_gamma, 'r')
                        
                        fluxint  = []
                        gammaint = []
                        
                        gt.hreadfilev4(file_flux,fluxint,param2)
                        gt.hreadfilev4(file_gamma,gammaint,param2)
                        
                        fluxint = fluxint[0]
                        gammaint = gammaint[0]
                        
                        file_flux.close()
                        file_gamma.close()
                        
                        int_flux = fluxint[0][0]
                        int_gamma = gammaint[0][0]
                        
                    except:
                        print "Data not found. Calculating."
                        # will estimate the fluxes considering propagation without regeneration or NC.
                        datapath = "../data/SunOscProbabilities/"
                        mu_inter = DM.DMNeuFluxDetNoInt(ch,dm,sig,param2,use_old_data = True, datapath = datapath)
                        
                        # for compatibility
                        int_gamma = 0.0
                        
                        # integrating
                        int_flux = integrate.quad(mu_inter,1.0*param.GeV,dm)[0]
                        
                        gammaint = [int_gamma]
                        fluxint  = [int_flux]
                        
                        #begin saving flux
                        filename_flux  = "DM_ANA_flux_"+str(dm/param.GeV)+"_GeV_xsection_"+str(sig/param.cm**2)+"_cm2_osc_"+param2.name+"_channel_"+ch+".dat"
                        filename_gamma = "DM_ANA_gamma_"+str(dm/param.GeV)+"_GeV_xsection_"+str(sig/param.cm**2)+"_cm2_osc_"+param2.name+"_channel_"+ch+".dat"
                        
                        file_flux  = open(savedatapath+filename_flux, 'w')
                        file_gamma = open(savedatapath+filename_gamma, 'w')
                        
                        gt.hwritefile(file_flux ,fluxint ,param2)
                        gt.hwritefile(file_gamma,gammaint,param2)            
                        #end saving flux
                        file_flux.close()
                        file_gamma.close()
    
                
                flux_lim.append(int_flux*param.km**2*param.year)
                gamma_lim.append(int_gamma*param.km**3*param.year)
                
            p2flim.append(flux_lim)
            p2glim.append(gamma_lim)        
               
    ################# MAKING PLOTS ######################        
    
    # Making new limits plot
    
    ### FLUX PLOT ###
    
    my_nu_flux_lim_soft = flim[0]
    my_nu_flux_lim_hard = flim[1]
    
    smy_nu_flux_lim_soft = sflim[0]
    smy_nu_flux_lim_hard = sflim[1]
    
    print smy_nu_flux_lim_soft
    print smy_nu_flux_lim_hard
    
    if param2 != None:
        my2_nu_flux_lim_soft = p2flim[0]
        my2_nu_flux_lim_hard = p2flim[1]
    
    # from Ref.: arXiv 0902.2460
    #DM_mass     = [250.0,500.0,1000.0]
    DM_mass     = [250.0,500.0,1000.0,3000.0,5000.0]
    
    #nu_flux_lim_hard = [8.8e2, 4.2e2, 3.6e2] # km^-2y^-1
    nu_flux_lim_hard = [8.8e2, 4.2e2, 3.6e2,3.3e2,3.6e2] # km^-2y^-1
    #nu_flux_lim_soft = [0.0  , 3.5e3, 1.3e3] # km^-2y^-1
    nu_flux_lim_soft = [0.0 , 3.5e3, 1.3e3,7.9e2,6.7e2] # km^-2y^-1
    
    fig = plt.figure()
    
    param.name = "3+3"
    sparam.name = "STD"
    
    plt.plot(DM_mass,nu_flux_lim_hard,'s',linestyle = "solid" ,label = "Icecube (hard)", color = "red")
    plt.plot(DM_mass,nu_flux_lim_soft,'s',linestyle = "dashed",label = "Icecube (soft)", color = "red")
    
    plt.plot(DM_mass,my_nu_flux_lim_hard,'o',linestyle = "solid" ,label = param.name+" Icecube (hard)", color = "green")
    plt.plot(DM_mass,my_nu_flux_lim_soft,'o',linestyle = "dashed",label = param.name+" Icecube (soft)", color = "green")
    
    print my_nu_flux_lim_soft
    print my_nu_flux_lim_hard
    
    plt.plot(DM_mass,smy_nu_flux_lim_hard,'o',linestyle = "solid" ,label = sparam.name+" Icecube (hard)", color = "blue")
    plt.plot(DM_mass,smy_nu_flux_lim_soft,'o',linestyle = "dashed",label = sparam.name+" Icecube (soft)", color = "blue")
    

    plt.loglog()
    
    plt.xlim(10.0,1.0e4)
    plt.ylim(50.0,1.0e5)
    
    plt.ylabel(r'$\mathrm{Muon\ flux\ from\ the\ Sun\ [km^{-2}y^{-1}]}$')
    plt.xlabel(r'$\mathrm{Neutralino\ mass\ [GeV]}$')
    
    plt.legend(loc = 'upper left')
    
    path = "../plots/"
    if use_MC_data :
        filename = "PlotNewICFluxLimitMC.eps"
    else: 
        filename = "PlotNewICFluxLimitAnalitic.eps"
    
    plt.savefig(path+filename)

    ### SIG PLOT ###
    
    # ratio : R = flux_2+3/flux_STD
    
    smy_nu_flux_lim_soft[0] =1.0
    
    R_hard = [my_nu_flux_lim_hard[i]/smy_nu_flux_lim_hard[i] for i in range(len(DM_mass))]
    R_soft = [my_nu_flux_lim_soft[i]/smy_nu_flux_lim_soft[i] for i in range(len(DM_mass))]
    
    if param2 != None:
        R2_hard = [my2_nu_flux_lim_hard[i]/smy_nu_flux_lim_hard[i] for i in range(len(DM_mass))]
        R2_soft = [my2_nu_flux_lim_soft[i]/smy_nu_flux_lim_soft[i] for i in range(len(DM_mass))]        
    
    #DMsigSD_soft   = [0.0,2.6e-38,2.2e-38]
    DMsigSD_soft   = [0.0,2.6e-38,2.2e-38,7.2e-38,1.5e-37]    
    #DMsigSD_hard   = [2.8e-40,3.0e-40,8.7e-40]
    DMsigSD_hard   = [2.8e-40,3.0e-40,8.7e-40,9.9e-39,3.6e-38]    
    
    # new DMsig
    R_soft[0] = 1.0
    my_DMsigSD_hard = [DMsigSD_hard[i]*R_hard[i]**-1 for i in range(len(DM_mass))]
    my_DMsigSD_soft = [DMsigSD_soft[i]*R_soft[i]**-1 for i in range(len(DM_mass))]
    
    if param2 != None:
        R2_soft[0] = 1.0
        my2_DMsigSD_hard = [DMsigSD_hard[i]*R2_hard[i]**-1 for i in range(len(DM_mass))]
        my2_DMsigSD_soft = [DMsigSD_soft[i]*R2_soft[i]**-1 for i in range(len(DM_mass))]        
    
    fig = plt.figure(figsize = (10,8))
    
    mpl.rcParams['axes.labelsize'] = "xx-large"
    mpl.rcParams['xtick.labelsize'] = "xx-large"
    mpl.rcParams['ytick.labelsize'] = "xx-large"
    mpl.rcParams['legend.fontsize'] = "small"
    
    mpl.rcParams['font.size'] = 16

    # CDMS limit
    DM_mass_CDMS = [9.9406783292,12.1694602882,14.4612794445,18.2381974079,23.9797671264,33.2631844688,46.6928782917,66.7251264804,95.351639623,137.890670608,193.562715148,274.964410758,392.930051792,558.174543467,792.911663423,1126.36614004,1600.05299448,2286.51011328,3228.8206912,4559.47384417,6400.31796953,10301.960886]
    DM_sig_CDMS  = [4.89012981667e-34,8.94244567855e-35,2.61502429561e-35,9.14460705548e-36,3.91050647843e-36,2.18678368583e-36,1.59912898498e-36,1.36748424071e-36,1.46234061987e-36,1.71005354324e-36,2.13844008852e-36,2.79641705226e-36,3.73951719026e-36,5.11373082554e-36,7.15103546736e-36,1e-35,1.43001240812e-35,1.99972775222e-35,2.73459618688e-35,3.99891108301e-35,5.5920727862e-35,8.94244567855e-35]
    plt.text(20.0,1.3e-36,u'CDMS (2008)', rotation=-32.0, color = "blue", fontsize = "large")
    
    # KIMS limit
    DM_mass_KIMS = [11.398495844,12.61174537,13.9541324973,15.8112573791,18.3470350855,22.0632380418,27.3333275932,35.5130192606,46.9715212036,63.2460768946,83.6528265901,111.96843479,148.979535199,199.407552045,265.321601349,350.929306576,461.405415313,613.922704562,807.191805189,1080.4178016,1420.54429515,1867.74606222,2485.12842698,3286.97094846,4373.47726647,5716.18061961,7651.04806427,10301.960886]
    DM_sig_KIMS  = [9.77892830637e-35,3.73951719026e-35,1.30769015567e-35,5.0006807121e-36,1.91228843285e-36,8.55143177036e-37,4.27629799149e-37,2.67414200583e-37,2.09116523131e-37,1.87001314859e-37,1.8286724512e-37,1.87001314859e-37,2.09116523131e-37,2.33847127116e-37,2.73459618688e-37,3.34404293379e-37,4.27629799149e-37,5.46844788602e-37,6.83835206662e-37,8.94244567855e-37,1.16939481815e-36,1.52920609179e-36,1.99972775222e-36,2.61502429561e-36,3.34404293379e-36,4.47183156237e-36,5.84777011194e-36,8.55143177036e-36]
    plt.text(1000.0,5.0e-37,u'KIMS (2007)', rotation=23.0, color = "red", fontsize = "large")
    
    # SUPER-K limit
    DM_mass_SUPERK = [18.209972014,21.5017707335,24.3555498225,27.5880909768,31.6227766017,37.3391948061,42.5467096332,50.8376617269,60.0275358048,68.8064047745,75.6599377056,78.4024974306,80.2858806092,82.2145063797,86.2118590704,88.2828385734,88.8083105364,90.4035671034,98.8201163355,113.272298009,126.042051687,141.925964953,161.719684948,178.886398103,200.237950626,239.257730776,279.174848197,321.908140042,382.361561747,446.153738108,517.508545714,636.9769417,743.248464437,846.905691618,976.541539256,1011.93971135]
    DM_sig_SUPERK  = [2.36541696568e-38,2.07707479418e-38,1.99030658719e-38,1.82532501751e-38,1.71132830416e-38,1.53603662455e-38,1.50458861493e-38,1.41118099726e-38,1.3827580497e-38,1.26828107778e-38,1.26942850868e-38,1.26985906285e-38,1.21564309444e-38,9.14312689313e-39,5.7715320858e-39,4.73887237089e-39,4.15486293154e-39,3.724013881e-39,3.64632747705e-39,3.49439992023e-39,3.27523905803e-39,3.20782090206e-39,3.0739903774e-39,3.01020431859e-39,4.99010823617e-39,5.10940292386e-39,5.11691668353e-39,5.59361646939e-39,6.11645457259e-39,6.68702868321e-39,7.80752709821e-39,9.52986816027e-39,1.13740612363e-38,1.35720686622e-38,1.65557719473e-38,1.73039132785e-38]
    plt.text(110.0,1.7e-39,u'SUPER-K (1996-2001)', rotation=0.0, color = "#FFA500", fontsize = "large")
    
    # PICASSO limit
    DM_mass_PICASSO = [4.72536,5.71749,6.74435,8.26474,10.79216,14.82722,20.63143,30.20445,47.72227,99.71682,248.92502,463.93097,2329.41406]
    DM_sig_PICASSO  = [3.54897E-35,1.41021E-35,6.60756E-36,3.53230E-36,2.22663E-36,1.60140E-36,1.23021E-36,1.15173E-36,1.27144E-36,1.95159E-36,4.16516E-36,7.53877E-36,3.32256E-35]
    plt.text(75.0,1.8e-36,u'PICASSO (2009)', rotation=22.0, color = "#b8860b", fontsize = "large")
    
    # XENON10 limit
    DM_mass_XENON10 = [8.22829,9.00321,10.15114,11.79400,13.70275,17.16041,22.14506,27.32009,34.21385,48.31023,92.08086,301.18085,1179.40027]
    DM_sig_XENON10  = [1.00000E-34,3.12071E-35,9.23671E-36,3.37859E-36,1.52724E-36,7.67463E-37,5.29832E-37,4.76608E-37,4.76608E-37,5.58633E-37,9.48444E-37,2.88251E-36,1.14149E-35]
    plt.text(45.0,3.0e-37,u'XENON10 (2008)', rotation=22.0, color = "magenta", fontsize = "large")
    
    # COUPP limit
    DM_mass_COUPP = [7.7161532813,8.0119016669,8.135964625,8.3213101601,8.6411585666,8.7734333512,9.1783630323,9.7459007869,9.9696630084,10.3503390196,10.7455505496,10.9907289255,11.9366119699,12.3018413306,12.7720139591,13.5598682997,14.6126112223,15.1711007722,16.3489344069,18.1503088137,20.9181054722,25.7708937217,30.8133252346,38.5064781618,49.5494574113,65.16861536,92.9480935779,124.923332872,162.987674449,223.914727698,356.513181554,435.144158599,563.42741112,785.399229049,1032.11016478,1619.04113494,2520.83251309,3839.23425575,5635.25802884,7459.40364987,9874.02928615]
    DM_sig_COUPP  = [1.03785718015e-36,8.63030145378e-37,7.25413688018e-37,6.57929047022e-37,5.29561865406e-37,4.69960147735e-37,3.78280401506e-37,2.97963567653e-37,2.5595990017e-37,2.2229410304e-37,1.93056288167e-37,1.73204972823e-37,1.22401315638e-37,1.04014314127e-37,8.93577469392e-38,7.35104654571e-38,5.85390594568e-38,5.02903711366e-38,4.00480531277e-38,3.0540138722e-38,2.06704738019e-38,1.36940699082e-38,1.00024449162e-38,7.88477245729e-39,6.92950040812e-39,6.78936693894e-39,7.49910577408e-39,8.74284723679e-39,1.04152440867e-38,1.32462711421e-38,2.00581890559e-38,2.41484182955e-38,3.07037586034e-38,4.21370795479e-38,5.35793719787e-38,8.29085240847e-38,1.31105044273e-37,1.98484446208e-37,2.90808432152e-37,3.86209530769e-37,5.12907416587e-37]
    plt.text(10.5,2.0e-38,u'COUPP (2011)', rotation=-56.0, color = "cyan", fontsize = "large")
    
    # SIMPLE limit
    DM_mass_SIMPLE = [5.1295317708,5.4471240212,5.6099146437,5.9101280574,6.3697469604,6.5626038562,6.9694536254,7.6299697108,8.3549894666,9.3588103613,10.8946136872,11.7552584534,14.1157981884,18.8760945376,25.0563307887,37.627429116,56.5292391067,79.8926160211,109.506423414,141.163345577,196.518835314,275.720663345,411.465635407,559.841578417,797.834060296,1077.23576584,1410.45391394,1790.90885197,2057.02463286]
    DM_sig_SIMPLE  = [1.02570691691e-36,6.88470147337e-37,4.72696175821e-37,2.92936426335e-37,1.67633239559e-37,1.2899203092e-37,8.85779986387e-38,5.55328855019e-38,3.72801865757e-38,2.25891676887e-38,1.41663060924e-38,1.14119716085e-38,7.9311885194e-39,5.77230402549e-39,4.60201715721e-39,4.35582027069e-39,4.67356205951e-39,5.49159719973e-39,6.60061620959e-39,8.02211437544e-39,1.02084673453e-38,1.34431322245e-38,1.93988965041e-38,2.55417673486e-38,3.6431572083e-38,4.79662016579e-38,6.24275313359e-38,8.12364592284e-38,9.11074283181e-38]
    plt.text(18.0,2.5e-39,u'SIMPLE (2011)', rotation=0.0, color = "green", fontsize = "large")
    
    
    # plotting other experiments limits
    
    plt.plot(DM_mass_CDMS,DM_sig_CDMS,      linestyle = "solid" ,label = "CDMS (2008)", color = "blue",lw = 4)
    plt.plot(DM_mass_KIMS,DM_sig_KIMS,      linestyle = "solid" ,label = "KIMS (2007)", color = "red",lw = 4)
    plt.plot(DM_mass_SUPERK,DM_sig_SUPERK,  linestyle = "solid" ,label = "SUPER-K (1996-2001)", color = "#FFA500",lw = 4)
    plt.plot(DM_mass_PICASSO,DM_sig_PICASSO,linestyle = "solid" ,label = "PICASSO (2009)", color = "#b8860b",lw = 4)
    plt.plot(DM_mass_XENON10,DM_sig_XENON10,linestyle = "solid" ,label = "XENON10 (2008)", color = "magenta",lw = 4)
    plt.plot(DM_mass_COUPP,DM_sig_COUPP,    linestyle = "solid" ,label = "COUPP (2011)", color = "cyan",lw = 4)
    plt.plot(DM_mass_SIMPLE,DM_sig_SIMPLE,  linestyle = "solid" ,label = "SIMPLE (2011)", color = "green",lw = 4)
    
    #icecube limits 
    plt.plot(DM_mass,DMsigSD_hard,'s',linestyle = "solid" ,label = r"STD Icecube ($\xi \bar{\xi} \to W^+ W^- $)", color = "k",lw = 4)
    plt.plot(DM_mass,my_DMsigSD_hard,'o',linestyle = "dashed" ,label = param.name+" Icecube ($\\xi \\bar{\\xi} \\to W^+ W^- $)", color = "k",lw = 4)
    if param2 != None :
        plt.plot(DM_mass,my2_DMsigSD_hard,'o',linestyle = "dotted" ,label = param2.name+" Icecube ($\\xi \\bar{\\xi} \\to W^+ W^- $)", color = "k",lw = 4)
    #plt.text(1100.0,0.6e-39,r"Icecube ($\chi \bar{\chi} \to W^+ W^- $)", rotation=46.0, color = "k", fontsize = "large")
    plt.text(150.0,1.5e-40,r"Icecube ($\chi \bar{\chi} \to W^+ W^- $)", rotation=0.0, color = "k", fontsize = "large")
    
    plt.plot(DM_mass,DMsigSD_soft,'s',linestyle = "solid",label = r"STD Icecube ($\xi \bar{\xi} \to b \bar{b} $)", color = "k")
    plt.plot(DM_mass,my_DMsigSD_soft,'o',linestyle = "dashed",label = param.name+" Icecube ($\\xi \\bar{\\xi} \\to b \\bar{b} $)", color = "k",lw = 4)
    if param2 != None :
        plt.plot(DM_mass,my2_DMsigSD_soft,'o',linestyle = "dotted",label = param2.name+" Icecube ($\\xi \\bar{\\xi} \\to b \\bar{b} $)", color = "k",lw = 4)
    #plt.text(2000.0,2.75e-38,r"Icecube ($\chi \bar{\chi} \to b \bar{b} $)", rotation=30.0, color = "k", fontsize = "large")
    plt.text(1000.0,2.75e-38,r"Icecube ($\chi \bar{\chi} \to b \bar{b} $)", rotation=0.0, color = "k", fontsize = "large")
    
    

    plt.loglog()
    
    plt.xlim(10.0,1.0e4)
    plt.ylim(1.0e-40,1.0e-35)
    
    plt.ylabel('$\mathrm{Dark\ Matter-proton\ SD\ cross\ section\ [cm^{2}]}$')
    plt.xlabel('$\mathrm{Dark\ Matter\ mass\ [GeV]}$')
    
    fig.subplots_adjust(bottom = 0.12, top = 0.95, left = 0.15, right = 0.95)
    
    #mpl.rcParams['legend.fontsize'] = "small"
    #plt.legend(loc = 'lower left',fancybox = True)
    
    path = "../plots/"
    if use_MC_data :
        filename = "PlotNewICXsectionLimitMC.eps"
    else: 
        filename = "PlotNewICXsectionLimitAnalitic.eps"
    
    
    plt.savefig(path+filename)    
    
def PlotIcecubeLimitCompareMCwithANA(param):
    """Plots the fluxes calculated using analytic calculation and monte carlo, and compares them.
    # param         : physics parameter
    """
    #DMm     = [250.0*param.GeV,500.0*param.GeV,1000.0*param.GeV]
    DMm     = [250.0*param.GeV,500.0*param.GeV,1000.0*param.GeV,3000.0*param.GeV,5000.0*param.GeV]
    DMsig_soft   = [0.0,2.5e-41*param.cm**2,1.8e-41*param.cm**2,5.3e-41*param.cm**2,1.1e-40*param.cm**2]
    DMsig_hard   = [3.70e-43*param.cm**2,2.9e-43*param.cm**2,7.2e-43*param.cm**2,7.4e-42*param.cm**2,2.6e-41*param.cm**2]

    savedatapath = "../data/DMConvolutedFlux/test/"
    
    #channels = ['WW','bb']
    channels = ['bb','WW']
    
    glim_MC = []
    flim_MC = []
    
    glim_ANA = []
    flim_ANA = []
    
    for ch in channels : 
        if ch == 'bb':
            DMsig = DMsig_soft
        elif ch == 'WW':
            DMsig = DMsig_hard
        else :
            print "Missing cross sections for this channel."
            quit()
            
        dat_pair = [[DMm[i],DMsig[i]] for i in range(len(DMm))]        
        
        E_nu_list_Full  = gt.LogSpaceEnergies(1.0,10000.0,binnum = 200)
        E_nu_hpl_Full   = gt.MidPoint(gt.LogSpaceEnergies(1.0,10000.0,binnum = 200))    
        
        ########################## FOR PARAM ###################################
        
        param.Refresh()
        
        gamma_lim_MC = []
        flux_lim_MC  = []
        
        gamma_lim_ANA = []
        flux_lim_ANA  = []        
        
        for dm, sig in dat_pair:
            ## USING MC DATA ##
            try :
                print "reading data"
                filename_flux  = "DM_flux_"+str(dm/param.GeV)+"_GeV_xsection_"+str(sig/param.cm**2)+"_cm2_osc_"+param.name+"_channel_"+ch+".dat"
                filename_gamma = "DM_gamma_"+str(dm/param.GeV)+"_GeV_xsection_"+str(sig/param.cm**2)+"_cm2_osc_"+param.name+"_channel_"+ch+".dat"
                
                file_flux  = open(savedatapath+filename_flux, 'r')
                file_gamma = open(savedatapath+filename_gamma, 'r')
                
                fluxint  = []
                gammaint = []
                
                gt.hreadfilev4(file_flux,fluxint,param)
                gt.hreadfilev4(file_gamma,gammaint,param)
                
                file_flux.close()
                file_gamma.close()
                               
                fluxint = fluxint[0]
                gammaint = gammaint[0]
                
                file_flux.close()
                file_gamma.close()
            except :
                print "generating data"
         
                E_nu_list  = filter(lambda x : x<dm/param.GeV, E_nu_list_Full) # [GeV]
                E_nu_hpl   = filter(lambda x : x<=dm/param.GeV, E_nu_hpl_Full) # [GeV]
         
                param.neutype = "neutrino"
                datapath = "../data/myMC/"+param.name+"/"+param.neutype+"/"
                
                inter_neu,inter_aneu = DM.DMFNeuFluxMCDetv2(ch,dm,sig,param,use_old_data = True, datapath = datapath)
                
                fluxint_neu  = [[E,float(inter_neu(E))*xs.signuNCC(E,0)*param.cm**2*(0.918*param.gr*param.cm**-3)/(939.27*param.MeV)*ice.MuonRange(ice.MuonEnergy(E*param.GeV,0,param),param)*param.meter] for E in E_nu_hpl]
                gammaint_neu = [[E,float(inter_neu(E))*xs.signuNCC(E,0)*param.cm**2*(0.918*param.gr*param.cm**-3)/(939.27*param.MeV)] for E in E_nu_hpl]                    
                
                param.neutype = "antineutrino"
                datapath = "../data/myMC/"+param.name+"/"+param.neutype+"/"
                
                ainter_neu,ainter_aneu = DM.DMFNeuFluxMCDetv2(ch,dm,sig,param,use_old_data = True, datapath = datapath)                    
                
                fluxint_aneu  = [[E,float(ainter_neu(E))*xs.signuNCC(E,1)*param.cm**2*(0.918*param.gr*param.cm**-3)/(939.27*param.MeV)*ice.MuonRange(ice.MuonEnergy(E*param.GeV,1,param),param)*param.meter] for E in E_nu_hpl]
                gammaint_aneu = [[E,float(ainter_neu(E))*xs.signuNCC(E,1)*param.cm**2*(0.918*param.gr*param.cm**-3)/(939.27*param.MeV)] for E in E_nu_hpl]
                
                fluxint  = [[E_nu_hpl[i],fluxint_neu[i][1]+fluxint_aneu[i][1]] for i in range(len(E_nu_hpl))]
                gammaint = [[E_nu_hpl[i],gammaint_neu[i][1]+gammaint_aneu[i][1]] for i in range(len(E_nu_hpl))]
                  
                #begin saving flux
                filename_flux  = "DM_flux_"+str(dm/param.GeV)+"_GeV_xsection_"+str(sig/param.cm**2)+"_cm2_osc_"+param.name+"_channel_"+ch+".dat"
                filename_gamma = "DM_gamma_"+str(dm/param.GeV)+"_GeV_xsection_"+str(sig/param.cm**2)+"_cm2_osc_"+param.name+"_channel_"+ch+".dat"
                
                file_flux  = open(savedatapath+filename_flux, 'w')
                file_gamma = open(savedatapath+filename_gamma, 'w')
                
                gt.hwritefile(file_flux ,fluxint ,param)
                gt.hwritefile(file_gamma,gammaint,param)            
                #end saving flux
                file_flux.close()
                file_gamma.close()
                                    
            pflux  = [f[1] for f in fluxint]
            pgamma = [f[1] for f in gammaint]
            
            int_flux_MC  = sum(pflux)
            int_gamma_MC = sum(pgamma)
            
            flux_lim_MC.append(int_flux_MC)
            gamma_lim_MC.append(int_gamma_MC)            
                
            ## USING ANALITIC CALCULATION ##
            try:
                print "reading data"
                filename_flux  = "DM_ANA_flux_"+str(dm/param.GeV)+"_GeV_xsection_"+str(sig/param.cm**2)+"_cm2_osc_"+param.name+"_channel_"+ch+".dat"
                filename_gamma = "DM_ANA_gamma_"+str(dm/param.GeV)+"_GeV_xsection_"+str(sig/param.cm**2)+"_cm2_osc_"+param.name+"_channel_"+ch+".dat"
                
                file_flux  = open(savedatapath+filename_flux, 'r')
                file_gamma = open(savedatapath+filename_gamma, 'r')
                
                fluxint  = []
                gammaint = []
                
                gt.hreadfilev4(file_flux,fluxint,param)
                gt.hreadfilev4(file_gamma,gammaint,param)
                
                file_flux.close()
                file_gamma.close()                
                
                fluxint = fluxint[0]
                gammaint = gammaint[0]
                
                file_flux.close()
                file_gamma.close()
                
                int_flux = fluxint[0][0]
                int_gamma = gammaint[0][0]
                
            except:
                print "Data not found. Calculating."
                # will estimate the fluxes considering propagation without regeneration or NC.
                datapath = "../data/SunOscProbabilities/test/"
                mu_inter = DM.DMNeuFluxDetNoInt(ch,dm,sig,param,onlyosc = False, datapath = datapath)
                
                # for compatibility
                int_gamma = 0.0
                
                # integrating
                int_flux = integrate.quad(mu_inter,1.0*param.GeV,dm)[0]
                
                gammaint = [int_gamma]
                fluxint  = [int_flux]
                
                #begin saving flux
                filename_flux  = "DM_ANA_flux_"+str(dm/param.GeV)+"_GeV_xsection_"+str(sig/param.cm**2)+"_cm2_osc_"+param.name+"_channel_"+ch+".dat"
                filename_gamma = "DM_ANA_gamma_"+str(dm/param.GeV)+"_GeV_xsection_"+str(sig/param.cm**2)+"_cm2_osc_"+param.name+"_channel_"+ch+".dat"
                
                file_flux  = open(savedatapath+filename_flux, 'w')
                file_gamma = open(savedatapath+filename_gamma, 'w')
                
                gt.hwritefile(file_flux ,fluxint ,param)
                gt.hwritefile(file_gamma,gammaint,param)            
                #end saving flux
                file_flux.close()
                file_gamma.close()

            
            flux_lim_ANA.append(int_flux*param.km**2*param.year)
            gamma_lim_ANA.append(int_gamma*param.km**3*param.year)
            
        # mc calculation
        flim_MC.append(flux_lim_MC)
        glim_MC.append(gamma_lim_MC)
        
        # analitic calculation
        flim_ANA.append(flux_lim_ANA)
        glim_ANA.append(gamma_lim_ANA)
               
    ################# MAKING PLOTS ######################        
    
    # Making new limits plot
    
    ### FLUX PLOT ###
    
    MC_nu_flux_lim_soft = map(lambda x : x*param.km**2*param.year,flim_MC[0])
    MC_nu_flux_lim_hard = map(lambda x : x*param.km**2*param.year,flim_MC[1])
    
    print MC_nu_flux_lim_hard

    ANA_nu_flux_lim_soft = flim_ANA[0]
    ANA_nu_flux_lim_hard = flim_ANA[1]
    
    # from Ref.: arXiv 0902.2460
    DM_mass     = [250.0,500.0,1000.0,3000.0,5000.0]
    
    #nu_flux_lim_hard = [8.8e2, 4.2e2, 3.6e2] # km^-2y^-1
    nu_flux_lim_hard = [8.8e2, 4.2e2, 3.6e2,3.3e2,3.6e2] # km^-2y^-1
    #nu_flux_lim_soft = [0.0  , 3.5e3, 1.3e3] # km^-2y^-1
    nu_flux_lim_soft = [0.0 , 3.5e3, 1.3e3,7.9e2,6.7e2] # km^-2y^-1
    
    fig = plt.figure(figsize = (10,10))
    
    mpl.rcParams['axes.labelsize'] = "x-large"
    mpl.rcParams['xtick.labelsize'] = "x-large"
    mpl.rcParams['ytick.labelsize'] = "x-large"
    mpl.rcParams['legend.fontsize'] = "small"
    
    mpl.rcParams['font.size'] = 22   
    
    plt.plot(DM_mass,nu_flux_lim_hard,'s',linestyle = "solid" ,label = "Icecube (hard)", color = "red")
    plt.plot(DM_mass,nu_flux_lim_soft,'s',linestyle = "dashed",label = "Icecube (soft)", color = "red")
    
    plt.plot(DM_mass,MC_nu_flux_lim_hard,'o',linestyle = "solid" ,label = param.name+"-MC Icecube (hard)", color = "green")
    plt.plot(DM_mass,MC_nu_flux_lim_soft,'o',linestyle = "dashed",label = param.name+"-MC Icecube (soft)", color = "green")
    
    plt.plot(DM_mass,ANA_nu_flux_lim_hard,'o',linestyle = "solid" ,label = param.name+"-ANA Icecube (hard)", color = "blue")
    plt.plot(DM_mass,ANA_nu_flux_lim_soft,'o',linestyle = "dashed",label = param.name+"-ANA Icecube (soft)", color = "blue")    

    plt.loglog()
    
    plt.xlim(10.0,1.0e4)
    plt.ylim(50.0,1.0e5)
    
    plt.ylabel('$\mathrm{Muon\ flux\ from\ the\ Sun\ [km^{-2}y^{-1}]}$')
    plt.xlabel('$\mathrm{Neutralino\ mass\ [GeV]}$')
    
    plt.legend(loc = 'upper left')
    
    path = "../plots/"
    
    filename = "PlotCompareICFluxLimitMCwithANA.eps"
    
    plt.savefig(path+filename)


    quit()
    ### SIG PLOT ###
    
    # ratio : R = flux_2+3/flux_STD
    
    ANA_nu_flux_lim_soft[0] =1.0
    MC_nu_flux_lim_soft[0] =1.0
    
    R_hard = [my_nu_flux_lim_hard[i]/smy_nu_flux_lim_hard[i] for i in range(len(DM_mass))]
    R_soft = [my_nu_flux_lim_soft[i]/smy_nu_flux_lim_soft[i] for i in range(len(DM_mass))]
    
    if param2 != None:
        R2_hard = [my2_nu_flux_lim_hard[i]/smy_nu_flux_lim_hard[i] for i in range(len(DM_mass))]
        R2_soft = [my2_nu_flux_lim_soft[i]/smy_nu_flux_lim_soft[i] for i in range(len(DM_mass))]        
    
    #DMsigSD_soft   = [0.0,2.6e-38,2.2e-38]
    DMsigSD_soft   = [0.0,2.6e-38,2.2e-38,7.2e-38,1.5e-37]    
    #DMsigSD_hard   = [2.8e-40,3.0e-40,8.7e-40]
    DMsigSD_hard   = [2.8e-40,3.0e-40,8.7e-40,9.9e-39,3.6e-38]    
    
    # new DMsig
    R_soft[0] = 1.0
    my_DMsigSD_hard = [DMsigSD_hard[i]*R_hard[i]**-1 for i in range(len(DM_mass))]
    my_DMsigSD_soft = [DMsigSD_soft[i]*R_soft[i]**-1 for i in range(len(DM_mass))]
    
    if param2 != None:
        R2_soft[0] = 1.0
        my2_DMsigSD_hard = [DMsigSD_hard[i]*R2_hard[i]**-1 for i in range(len(DM_mass))]
        my2_DMsigSD_soft = [DMsigSD_soft[i]*R2_soft[i]**-1 for i in range(len(DM_mass))]        
    
    fig = plt.figure()
    
    #icecube limits 
    plt.plot(DM_mass,DMsigSD_hard,'s',linestyle = "solid" ,label = r"STD Icecube ($\xi \bar{\xi} \to W^+ W^- $)", color = "k")
    plt.plot(DM_mass,my_DMsigSD_hard,'o',linestyle = "dashed" ,label = param.name+" Icecube ($\\xi \\bar{\\xi} \\to W^+ W^- $)", color = "k")
    if param2 != None :
        plt.plot(DM_mass,my2_DMsigSD_hard,'o',linestyle = "dotted" ,label = param2.name+" Icecube ($\\xi \\bar{\\xi} \\to W^+ W^- $)", color = "k")
    plt.text(1100.0,0.6e-39,r"Icecube ($\chi \bar{\chi} \to W^+ W^- $)", rotation=46.0, color = "k", fontsize = "large")
    
    plt.plot(DM_mass,DMsigSD_soft,'s',linestyle = "solid",label = r"STD Icecube ($\xi \bar{\xi} \to b \bar{b} $)", color = "k")
    plt.plot(DM_mass,my_DMsigSD_soft,'o',linestyle = "dashed",label = param.name+" Icecube ($\\xi \\bar{\\xi} \\to b \\bar{b} $)", color = "k")
    if param2 != None :
        plt.plot(DM_mass,my2_DMsigSD_soft,'o',linestyle = "dotted",label = param2.name+" Icecube ($\\xi \\bar{\\xi} \\to b \\bar{b} $)", color = "k")
    plt.text(2000.0,2.75e-38,r"Icecube ($\chi \bar{\chi} \to b \bar{b} $)", rotation=30.0, color = "k", fontsize = "large")
    
    # CDMS limit
    DM_mass_CDMS = [9.9406783292,12.1694602882,14.4612794445,18.2381974079,23.9797671264,33.2631844688,46.6928782917,66.7251264804,95.351639623,137.890670608,193.562715148,274.964410758,392.930051792,558.174543467,792.911663423,1126.36614004,1600.05299448,2286.51011328,3228.8206912,4559.47384417,6400.31796953,10301.960886]
    DM_sig_CDMS  = [4.89012981667e-34,8.94244567855e-35,2.61502429561e-35,9.14460705548e-36,3.91050647843e-36,2.18678368583e-36,1.59912898498e-36,1.36748424071e-36,1.46234061987e-36,1.71005354324e-36,2.13844008852e-36,2.79641705226e-36,3.73951719026e-36,5.11373082554e-36,7.15103546736e-36,1e-35,1.43001240812e-35,1.99972775222e-35,2.73459618688e-35,3.99891108301e-35,5.5920727862e-35,8.94244567855e-35]
    plt.text(21.0,1.75e-36,u'CDMS (2008)', rotation=-40.0, color = "blue", fontsize = "large")
    
    # KIMS limit
    DM_mass_KIMS = [11.398495844,12.61174537,13.9541324973,15.8112573791,18.3470350855,22.0632380418,27.3333275932,35.5130192606,46.9715212036,63.2460768946,83.6528265901,111.96843479,148.979535199,199.407552045,265.321601349,350.929306576,461.405415313,613.922704562,807.191805189,1080.4178016,1420.54429515,1867.74606222,2485.12842698,3286.97094846,4373.47726647,5716.18061961,7651.04806427,10301.960886]
    DM_sig_KIMS  = [9.77892830637e-35,3.73951719026e-35,1.30769015567e-35,5.0006807121e-36,1.91228843285e-36,8.55143177036e-37,4.27629799149e-37,2.67414200583e-37,2.09116523131e-37,1.87001314859e-37,1.8286724512e-37,1.87001314859e-37,2.09116523131e-37,2.33847127116e-37,2.73459618688e-37,3.34404293379e-37,4.27629799149e-37,5.46844788602e-37,6.83835206662e-37,8.94244567855e-37,1.16939481815e-36,1.52920609179e-36,1.99972775222e-36,2.61502429561e-36,3.34404293379e-36,4.47183156237e-36,5.84777011194e-36,8.55143177036e-36]
    plt.text(1000.0,5.5e-37,u'KIMS (2007)', rotation=23.0, color = "red", fontsize = "large")
    
    # SUPER-K limit
    DM_mass_SUPERK = [18.209972014,21.5017707335,24.3555498225,27.5880909768,31.6227766017,37.3391948061,42.5467096332,50.8376617269,60.0275358048,68.8064047745,75.6599377056,78.4024974306,80.2858806092,82.2145063797,86.2118590704,88.2828385734,88.8083105364,90.4035671034,98.8201163355,113.272298009,126.042051687,141.925964953,161.719684948,178.886398103,200.237950626,239.257730776,279.174848197,321.908140042,382.361561747,446.153738108,517.508545714,636.9769417,743.248464437,846.905691618,976.541539256,1011.93971135]
    DM_sig_SUPERK  = [2.36541696568e-38,2.07707479418e-38,1.99030658719e-38,1.82532501751e-38,1.71132830416e-38,1.53603662455e-38,1.50458861493e-38,1.41118099726e-38,1.3827580497e-38,1.26828107778e-38,1.26942850868e-38,1.26985906285e-38,1.21564309444e-38,9.14312689313e-39,5.7715320858e-39,4.73887237089e-39,4.15486293154e-39,3.724013881e-39,3.64632747705e-39,3.49439992023e-39,3.27523905803e-39,3.20782090206e-39,3.0739903774e-39,3.01020431859e-39,4.99010823617e-39,5.10940292386e-39,5.11691668353e-39,5.59361646939e-39,6.11645457259e-39,6.68702868321e-39,7.80752709821e-39,9.52986816027e-39,1.13740612363e-38,1.35720686622e-38,1.65557719473e-38,1.73039132785e-38]
    plt.text(110.0,2.0e-39,u'SUPER-K (1996-2001)', rotation=0.0, color = "#FFA500", fontsize = "large")
    
    # PICASSO limit
    DM_mass_PICASSO = [4.72536,5.71749,6.74435,8.26474,10.79216,14.82722,20.63143,30.20445,47.72227,99.71682,248.92502,463.93097,2329.41406]
    DM_sig_PICASSO  = [3.54897E-35,1.41021E-35,6.60756E-36,3.53230E-36,2.22663E-36,1.60140E-36,1.23021E-36,1.15173E-36,1.27144E-36,1.95159E-36,4.16516E-36,7.53877E-36,3.32256E-35]
    plt.text(120.0,2.5e-36,u'PICASSO (2009)', rotation=22.0, color = "#b8860b", fontsize = "large")
    
    # XENON10 limit
    DM_mass_XENON10 = [8.22829,9.00321,10.15114,11.79400,13.70275,17.16041,22.14506,27.32009,34.21385,48.31023,92.08086,301.18085,1179.40027]
    DM_sig_XENON10  = [1.00000E-34,3.12071E-35,9.23671E-36,3.37859E-36,1.52724E-36,7.67463E-37,5.29832E-37,4.76608E-37,4.76608E-37,5.58633E-37,9.48444E-37,2.88251E-36,1.14149E-35]
    plt.text(45.0,3.5e-37,u'XENON10 (2008)', rotation=22.0, color = "magenta", fontsize = "large")
    
    # COUPP limit
    DM_mass_COUPP = [7.7161532813,8.0119016669,8.135964625,8.3213101601,8.6411585666,8.7734333512,9.1783630323,9.7459007869,9.9696630084,10.3503390196,10.7455505496,10.9907289255,11.9366119699,12.3018413306,12.7720139591,13.5598682997,14.6126112223,15.1711007722,16.3489344069,18.1503088137,20.9181054722,25.7708937217,30.8133252346,38.5064781618,49.5494574113,65.16861536,92.9480935779,124.923332872,162.987674449,223.914727698,356.513181554,435.144158599,563.42741112,785.399229049,1032.11016478,1619.04113494,2520.83251309,3839.23425575,5635.25802884,7459.40364987,9874.02928615]
    DM_sig_COUPP  = [1.03785718015e-36,8.63030145378e-37,7.25413688018e-37,6.57929047022e-37,5.29561865406e-37,4.69960147735e-37,3.78280401506e-37,2.97963567653e-37,2.5595990017e-37,2.2229410304e-37,1.93056288167e-37,1.73204972823e-37,1.22401315638e-37,1.04014314127e-37,8.93577469392e-38,7.35104654571e-38,5.85390594568e-38,5.02903711366e-38,4.00480531277e-38,3.0540138722e-38,2.06704738019e-38,1.36940699082e-38,1.00024449162e-38,7.88477245729e-39,6.92950040812e-39,6.78936693894e-39,7.49910577408e-39,8.74284723679e-39,1.04152440867e-38,1.32462711421e-38,2.00581890559e-38,2.41484182955e-38,3.07037586034e-38,4.21370795479e-38,5.35793719787e-38,8.29085240847e-38,1.31105044273e-37,1.98484446208e-37,2.90808432152e-37,3.86209530769e-37,5.12907416587e-37]
    plt.text(11.0,3.0e-38,u'COUPP (2011)', rotation=-56.0, color = "cyan", fontsize = "large")
    
    # SIMPLE limit
    DM_mass_SIMPLE = [5.1295317708,5.4471240212,5.6099146437,5.9101280574,6.3697469604,6.5626038562,6.9694536254,7.6299697108,8.3549894666,9.3588103613,10.8946136872,11.7552584534,14.1157981884,18.8760945376,25.0563307887,37.627429116,56.5292391067,79.8926160211,109.506423414,141.163345577,196.518835314,275.720663345,411.465635407,559.841578417,797.834060296,1077.23576584,1410.45391394,1790.90885197,2057.02463286]
    DM_sig_SIMPLE  = [1.02570691691e-36,6.88470147337e-37,4.72696175821e-37,2.92936426335e-37,1.67633239559e-37,1.2899203092e-37,8.85779986387e-38,5.55328855019e-38,3.72801865757e-38,2.25891676887e-38,1.41663060924e-38,1.14119716085e-38,7.9311885194e-39,5.77230402549e-39,4.60201715721e-39,4.35582027069e-39,4.67356205951e-39,5.49159719973e-39,6.60061620959e-39,8.02211437544e-39,1.02084673453e-38,1.34431322245e-38,1.93988965041e-38,2.55417673486e-38,3.6431572083e-38,4.79662016579e-38,6.24275313359e-38,8.12364592284e-38,9.11074283181e-38]
    plt.text(22.0,3.0e-39,u'SIMPLE (2011)', rotation=0.0, color = "green", fontsize = "large")
    
    
    # plotting other experiments limits
    
    plt.plot(DM_mass_CDMS,DM_sig_CDMS,      linestyle = "solid" ,label = "CDMS (2008)", color = "blue")
    plt.plot(DM_mass_KIMS,DM_sig_KIMS,      linestyle = "solid" ,label = "KIMS (2007)", color = "red")
    plt.plot(DM_mass_SUPERK,DM_sig_SUPERK,  linestyle = "solid" ,label = "SUPER-K (1996-2001)", color = "#FFA500")
    plt.plot(DM_mass_PICASSO,DM_sig_PICASSO,linestyle = "solid" ,label = "PICASSO (2009)", color = "#b8860b")
    plt.plot(DM_mass_XENON10,DM_sig_XENON10,linestyle = "solid" ,label = "XENON10 (2008)", color = "magenta")
    plt.plot(DM_mass_COUPP,DM_sig_COUPP,    linestyle = "solid" ,label = "COUPP (2011)", color = "cyan")
    plt.plot(DM_mass_SIMPLE,DM_sig_SIMPLE,  linestyle = "solid" ,label = "SIMPLE (2011)", color = "green")

    plt.loglog()
    
    plt.xlim(10.0,1.0e4)
    plt.ylim(1.0e-40,1.0e-35)
    
    plt.ylabel('$\mathrm{Neutralino-proton\ SD\ cross-section\ [cm^{2}]}$')
    plt.xlabel('$\mathrm{Dark\ Matter\ mass\ [GeV]}$')
    
    #mpl.rcParams['legend.fontsize'] = "small"
    #plt.legend(loc = 'lower left',fancybox = True)
    
    path = "../plots/"
    filename = "PlotCompareNewICXsectionLimit.eps"
    
    plt.savefig(path+filename)        
    
def PlotAnnihilationRatesLimit(param):
    """Plots new Icecube limit considering 2+3 case.
    # param         : physics parameter
    """
    DMm     = [250.0*param.GeV,500.0*param.GeV,1000.0*param.GeV]
    # from paper
    #DMann_soft   = [0.0,1.4e23/param.sec,3.0e22/param.sec]
    #DMann_hard   = [6.0e21/param.sec,1.6e21/param.sec,1.2e21/sec]
    
    DMann_soft   = [0.0,1.4e23,3.0e22]
    DMann_hard   = [6.0e21,1.6e21,1.2e21]
    
    DMsig_soft   = [0.0,2.5e-41*param.cm**2,1.8e-41*param.cm**2]
    DMsig_hard   = [3.70e-43*param.cm**2,2.9e-43*param.cm**2,7.2e-43*param.cm**2]
    
    dat_pair_soft = [[DMm[i],DMsig_soft[i]] for i in range(len(DMm))]
    dat_pair_hard = [[DMm[i],DMsig_hard[i]] for i in range(len(DMm))]
    
    my_DMann_hard = [np.sum(DM.DMSunAnnihilationRate(dm,sig,param)*param.sec) for dm,sig in dat_pair_hard]
    my_DMann_soft = [np.sum(DM.DMSunAnnihilationRate(dm,sig,param)*param.sec) for dm,sig in dat_pair_soft]
    
    for i,a in enumerate(DMann_hard):
        print a/my_DMann_hard[i]
        print DMann_soft[i]/my_DMann_soft[i]
    
    # my limits
    
    fig = plt.figure()
    
    DM_mass     = [250.0,500.0,1000.0]
    
    plt.plot(DM_mass,DMann_hard,'s',linestyle = "solid" ,label = "Icecube (hard)", color = "red")
    plt.plot(DM_mass,DMann_soft,'s',linestyle = "dashed",label = "Icecube (soft)", color = "red")
    
    plt.plot(DM_mass,my_DMann_hard,'o',linestyle = "solid" ,label = "My Rate (hard)", color = "blue")
    plt.plot(DM_mass,my_DMann_soft,'o',linestyle = "dashed",label = "My Rate (soft)", color = "blue")    

    plt.loglog()
    
    plt.xlim(10.0,1.0e4)
    #plt.ylim(50.0,1.0e5)
    
    plt.ylabel('$\mathrm{Annihilation\ Rate [s^{-1}]}$')
    plt.xlabel('$\mathrm{Neutralino\ mass\ [GeV]}$')
    
    plt.legend(loc = 'upper left')
    
    path = "../plots/"
    filename = "PlotAnnihilationRatesLimit.png"
    
    plt.savefig(path+filename)
    
def PlotParameterScan(datapath):
    """Using calculated data from CalculateParameaterScan will make the contour plot.
    
    @type   datapath:   string
    @param  datapath:   path to zip files.
    
    @rtype          :   plot
    @return         :   generates the countour plot
    """
    
    ratio_flux_soft = []
    ratio_flux_hard = []
    
    #begin reading data #
    filename_soft = "3+3_soft_ratio.dat"
    file = open(datapath+filename_soft,'r')
    ratio_flux_soft = gt.quickread(file)
    file.close()
    filename_hard = "3+3_hard_ratio.dat"
    file = open(datapath+filename_hard,'r')
    ratio_flux_hard = gt.quickread(file)
    file.close()
    #end reading data #
    
    # generate contour plot #
    fig = plt.figure()
    ax = plt.subplot(111)    
    
    R_soft_plane = sorted(ratio_flux_soft,key = lambda x : x[0])
    R_hard_plane = sorted(ratio_flux_hard,key = lambda x : x[0])
    
    R_soft_plane = sorted(ratio_flux_soft,key = lambda x : x[1])
    R_hard_plane = sorted(ratio_flux_hard,key = lambda x : x[1])
    
    th      = map(lambda x : x[0],R_soft_plane)
    #redefine th as sin^2 2th
    th      = map(lambda x : np.sin(2.0*x)**2,th)
    dmsq    = map(lambda x : x[1],R_soft_plane)
    
    Rsoft = map(lambda x : x[2],R_soft_plane)
    Rhard = map(lambda x : x[2],R_hard_plane)
    
    minth = np.min(th)
    maxth = np.max(th)
    mindm = np.min(dmsq)
    maxdm = np.max(dmsq)
    
    stepth = (maxth - minth)/30.0
    stepdm = (maxdm - mindm)/30.0
    
    nth = np.arange(minth,maxth,stepth)
    ndm = np.arange(mindm,maxdm,stepdm)
    
    nRsoft = pl.griddata(th,dmsq,Rsoft,nth,ndm,interp = 'linear')
    nRhard = pl.griddata(th,dmsq,Rhard,nth,ndm,interp = 'linear')
    
    #nRsoft = interpolate.griddata(th,dmsq,Rsoft,nth,ndm,method = 'linear')
    #nRhard = interpolate.griddata(th,dmsq,Rhard,nth,ndm,method = 'linear')
    
    contours_soft = np.array([1.25,1.5,1.75,2.0])
    contours_hard = np.array([1.25,1.5,1.75,2.0])
    #contours_hard = np.array([5.0,7.0,8.0,8.5])
    
    CS_soft = plt.contour(nth,ndm,nRsoft,contours_soft, colors = "blue",linestyles=['solid','dashed','dashdot','dotted'], label = r"$W^+W^-$")
    CS_hard = plt.contour(nth,ndm,nRhard,contours_hard, colors = "red",linestyles=['solid','dashed','dashdot','dotted'], label = r"$b\bar{b}$")
    
    plt.clabel(CS_soft, inline=1, fontsize=8, fmt='%1.2f')
    plt.clabel(CS_hard, inline=1, fontsize=8, fmt='%1.2f')
    
    plt.xlabel(r"$\sin^2(2\theta_s)$")
    plt.ylabel(r"$\Delta m^2_s \mathrm{[eV^2]}$")
    
    #extra labels
    
    # extra legend
    box1t = osb.TextArea(r"$b\bar{b}$", textprops=dict(color="b"))
    box1d = osb.DrawingArea(60, 20, 0, 0)
    el1 = ptc.Ellipse((10, 10), width=5, height=5, angle=0, fc="b", edgecolor = 'none')
    box1d.add_artist(el1)
    
    
    box2t = osb.TextArea(r"$W^+W^-$", textprops=dict(color="r"))
    #box2t = osb.TextArea("2+3", textprops=dict(color="k"))
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
                                 child=box1, pad=0.25,
                                 frameon=False,
                                 #bbox_to_anchor=(0., 1.02),
                                 #bbox_transform=ax.transAxes,
                                 borderpad=0.,
                                 )
    
    anchored_box2 = osb.AnchoredOffsetbox(loc=9,
                                 child=box2, pad=1.25,
                                 frameon=False,
                                 #bbox_to_anchor=(0., 1.02),
                                 #bbox_transform=ax.transAxes,
                                 borderpad=0.,
                                 )    
    
    ax.add_artist(anchored_box1)
    ax.add_artist(anchored_box2)
    
    plt.loglog()
    
    plt.xlim(1.0e-3,0.2)
    #plt.ylim()
    
    filename = "countour.eps"
    plt.savefig(datapath + filename)
    
    
#===============================================================================
# Testing
#===============================================================================
    
if __name__ == '__main__':    
	pass
