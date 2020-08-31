"""
Author  : C.A. Arguelles
Date    : 10/MAY/2011

This script models the NUMI beam.

Log:
- Modified on 23/FEB/2012 by C.Arguelles
    + Review the code to make it compatible with the neutrino
    commons library.
"""

# python standard models
import numpy as np
import scipy.interpolate as interpolate

#===============================================================================
# NUMI beam
#===============================================================================

def NUMIflux_binned(E,neu):
  """ Binned NUMI flux measured at the near detector.
  
  # REFERENCE ARXIV 0910.2201v2 pag. 18
  @type   E 	:	float
  @param  E	:	neutrino energy [GeV]
  @type   neu	:	integer
  @param  neu	:	neu : 0 : neutrino, neu : 1 : antineutrino

  @rtype	:	float
  @return	:	flux [Particle/GeV/cm^2/POT]
  """
  if neu == 0 :
    if E <=0 :
      flux = 0.0
    elif E > 0.0 and E <= 1.0 :
      flux = 1.00e3 ## FROM PRESENTATION RESCALED
    elif E > 1.0 and E <= 2.0 :
      flux = 2.2e4 ## FROM PRESENTATION RESCALED
    elif E > 2.0 and E <= 3.0 :
      flux = 5.57e4 ## FROM PRESENTATION RESCALED
    elif E>3.0 and E<=4.0 :
      flux = 8.05e4 ## FROM PRESENTATION RESCALED
    elif E>4.0 and E<=5.0 :
      flux = 3.06e4
    elif E>5.0 and E<=7.0 :
      flux = 9.07e3
    elif E>7.0 and E<=9.0 :
      flux = 5.18e3
    elif E>9.0 and E<=12.0 :
      flux = 3.21e3
    elif E>12.0 and E<=15.0 :
      flux = 1.94e3
    elif E>15.0 and E<=18.0 :
      flux = 1.09e3
    elif E>18.0 and E<=22.0 :
      flux = 629.0
    elif E>22.0 and E<=26.0 :
      flux = 348.0
    elif E>26.0 and E<=30.0 :
      flux = 200.0
    elif E>30.0 and E<=36.0 :
      flux = 119.0
    elif E>36.0 and E<=42.0 :
      flux = 72.2
    elif E>42.0 and E<=50.0 :
      flux = 51.6
    elif E>50.0:
      flux = 0.0
  elif neu == 1:
    if E <=0 :
      flux = 0.0
    elif E > 0.0 and E <= 0.5 :
      flux = 0.0 ## FROM PRESENTATION RESCALED
    elif E > 0.5 and E <= 1.0 :
      flux = 0.837209302326e3 ## FROM PRESENTATION RESCALED
    elif E > 1.0 and E <= 1.5 :
      flux = 3.13953488372e3 ## FROM PRESENTATION RESCALED
    elif E > 1.5 and E <= 2.0 :
      flux = 6.54069767442e3 ## FROM PRESENTATION RESCALED
    elif E > 2.0 and E <= 2.5 :
      flux = 9.94186046512e3 ## FROM PRESENTATION RESCALED
    elif E > 2.5 and E <= 3.0 :
      flux = 13.3953488372e3 ## FROM PRESENTATION RESCALED
    elif E > 3.0 and E <= 3.5 :
      flux = 15.488372093e3 ## FROM PRESENTATION RESCALED
    elif E > 3.5 and E <= 4.0 :
      flux = 15.0174418605e3 ## FROM PRESENTATION RESCALED
    elif E > 4.0 and E <= 4.5 :
      flux = 10.9360465116e3 ## FROM PRESENTATION RESCALED
    elif E > 4.5 and E <= 5.0 :
      flux = 6.43604651163e3 ## FROM PRESENTATION RESCALED
    elif E>5.0 and E<=7.0 :
      flux = 2.8e3
    elif E>7.0 and E<=9.0 :
      flux = 2.32e3
    elif E>9.0 and E<=12.0 :
      flux = 1.32e3
    elif E>12.0 and E<=15.0 :
      flux = 6.89e2
    elif E>15.0 and E<=18.0 :
      flux = 3.79e2
    elif E>18.0 and E<=22.0 :
      flux = 190.0
    elif E>22.0 and E<=26.0 :
      flux = 86.3
    elif E>26.0 and E<=30.0 :
      flux = 40.1
    elif E>30.0 and E<=36.0 :
      flux = 1.9
    elif E>36.0 and E<=42.0 :
      flux = 0.9
    elif E>42.0 and E<=50.0 :
      flux = 0.5
    elif E>50.0:
      flux = 0.0
  return flux*1.0e-9*1.0e-4 # to conv. to m -> cm and 1.0POT

def NUMIflux(E,neu):
  """ Interpolated NUMI flux measured at the near detector.
  
  # REFERENCE ARXIV 0910.2201v2 pag. 18
  @type   E 	:	float
  @param  E	:	neutrino energy [GeV]
  @type   neu	:	integer
  @param  neu	:	neu : 0 : neutrino, neu : 1 : antineutrino

  @rtype	:	float
  @return	:	flux [Particle/GeV/cm^2/POT]
  """  
  Ev = [0.5,1.0,1.5,2.5,3.5,4.5,6.0,8.0,10.5,13.5,16.5,20.0,24.0,28.0,33.0,39.0,46.0,50.0]
  Eav = [1.5,2.5,3.5,4.5,6.0,8.0,10.5,13.5,16.5,20.0,24.0,28.0,33.0,39.0,46.0,50.0]
  fluxneu = [1.0e3,1.0e4,2.2e4,5.57e4,8.05e4,3.06e4,9.07e3,5.18e3,3.21e3,1.94e3,1.09e3,629.0,348.0,200.0,119.0,72.2,51.6,51.6]
  fluxaneu = [0.0,0.0,0.0,2.8e3,2.32e3,1.32e3,6.89e2,3.79e2,190.0,86.3,40.1,1.9,0.9,0.5,0.5]
  if neu == 0 : 
    #inter = interpolate.interp1d(Ev,fluxneu)
    inter = interpolate.InterpolatedUnivariateSpline(Ev,fluxneu)
    if E<0.5 : 
      return 0.0
    elif E>46.0 :
      return 0.0
    else :
      return inter(E)*1.0e-9*1.0e-4# to conv. to m -> cm and 1.0POT
  elif neu == 1 : 
    inter = interpolate.interp1d(Ev,fluxaneu)
    if E<2.5 : 
      return 0.0
    elif E>46.0 :
      return 0.0
    else :
      return inter(E)*1.0e-9*1.0e-4# to conv. to m -> cm and 1.0POT
      
#===============================================================================
# CORRECTION TO FLUX SHAPE DUE TO PIONS DECAY ANGLES
#===============================================================================      
      
def FNratio_binned(E):
  """ Ratio between neutrino flux at the Far Detector and at Near Detector.
  
  Ref. MINOS Document 5694-v1. Fig. 4.3
  Tesis : Zarkov Pavlovich
  URL : http://minos-docdb.fnal.gov/cgi-bin/ShowDocument?docid=5694
  
  @type   E 	:	float
  @param  E	:	neutrino energy [GeV]

  @rtype	:	float
  @return	:	F/N [dimensionless]
  """  
  if E <=0 :
    fnr = 0.0
  elif E > 0.0 and E <= 0.5 :
    fnr = 1.56980937443
  elif E > 0.5 and E <= 1.0 :
    fnr = 1.68239230683
  elif E > 1.0 and E <= 1.5 :
    fnr = 1.59944451298
  elif E> 1.5 and E<=2.0 :
    fnr = 1.4945939881
  elif E>2.0 and E<=2.5 :
    fnr = 1.35689923574
  elif E>2.5 and E<=3.0 :
    fnr = 1.23797378825
  elif E>3.0 and E<=3.5 :
    fnr = 1.17535790024
  elif E>3.5 and E<=4.0 :
    fnr = 1.23006538251
  elif E>4.0 and E<=4.5 :
    fnr = 1.41616622315
  elif E>4.5 and E<=5.0 :
    fnr = 1.55534133432
  elif E>5.0 and E<=5.5 :
    fnr = 1.57093938166
  elif E>5.5 and E<=6.0 :
    fnr = 1.46452790065
  elif E>6.0 and E<=6.5 :
    fnr = 1.3737572329
  elif E>6.5 and E<=7.0 :
    fnr = 1.29550053164
  elif E>7.0 and E<=7.5 :
    fnr = 1.27981366277
  elif E>7.5 and E<=8.0 :
    fnr = 1.24848762549
  elif E>8.0 and E<=9.0 :
    fnr = 1.26719442633
  elif E>9.0 and E<=10.0 :
    fnr = 1.24833629993
  elif E>10.0 and E<=11.0 :
    fnr = 1.26389158134
  elif E>11.0 and E<=12.0 :
    fnr = 1.34670943259
  elif E>12.0 and E<=14.0 :
    fnr = 1.41853479724
  elif E>14.0 and E<=16.0 :
    fnr = 1.56070845743
  elif E>16.0 and E<=18.0 :
    fnr = 1.55427712082
  elif E>18.0:
    fnr = 1.55427712082
  return (fnr*1.0e-6)      
      
def FNratio(E):
  """ Ratio between neutrino flux at the Far Detector and at Near Detector.
  
  Ref. MINOS Document 5694-v1. Fig. 4.3
  Tesis : Zarkov Pavlovich
  URL : http://minos-docdb.fnal.gov/cgi-bin/ShowDocument?docid=5694
  
  @type   E 	:	float
  @param  E	:	neutrino energy [GeV]

  @rtype	:	float
  @return	:	F/N [dimensionless]
  """
  Enu = [0.0,0.25,0.75,1.25,1.75,2.25,2.75,3.25,3.75,4.25,4.75,5.25,5.75,6.25,6.75,7.25,7.75,8.5,9.5,10.5,11.5,13.0,15.0,17.0]
  FNr = [1.569,1.569,1.682,1.599,1.494,1.356,1.237,1.175,1.230,1.416,1.555,1.571,1.465,1.374,1.296,1.280,1.248,1.267,1.248,1.264,1.347,1.418,1.561,1.554]
  
  if E>=17.0 : 
    return 1.554*1.0e-6
  else:
    inter = interpolate.interp1d(Enu,FNr)
    return inter(E)*1.0e-6

      
if __name__ == '__main__':
  import matplotlib.pyplot as plt
  import numpy as np
  
  import neuosc as no
  import physicsconstants as PC
  import body as bd
  E = np.arange(0.1,10.0,0.05)
  neutype = 0
  Flux = [NUMIflux(EE,neutype) for EE in E]
  
  # osc parameters
  ineu = 0
  fneu = 0
  
  pc = PC.PhysicsConstants()
  
  pc.numneu = 5
  
  # set std ang to zero
  pc.th12 = 0.0
  pc.th23 = 0.0
  pc.th13 = 0.0
    
  pc.th14 = 0.129599 # higher increases resonances amplitud # if <0.1 no effect
  pc.th24 = 0.171848 # no effect on pee
  pc.th15 = 0.138442
  pc.th25 = 0.149991
  pc.th34 = 0.15#0.05 # no effect on pee
  pc.dm41sq = -0.47 # higher moves resonante to right / proportional
  pc.dm51sq = -0.87
  pc.delta1 = 0.0
  
  pc.Refresh()
  
  L = 840.0*pc.meter
  
  body  = bd.Earth()
  track = body.track(L,0,L)
  
  fM2 = no.flavorM2(pc)
  
  OFlux = [NUMIflux(EE,neutype)*no.OscProbConstantDensity(ineu,fneu,EE,track,body,fM2,pc) for EE in E]
  
  plt.plot(E,Flux)
  plt.plot(E,OFlux)
  plt.ylabel(r'$\phi [\mathrm{GeV}^{-1}\mathrm{cm}^{-2}\mathrm{POT}^{-1}]$')
  plt.xlabel(r'$E [\mathrm{GeV}]$')
  
  plt.title(r'NUMI flux in a 2+3 scheme without STD-osc')
  plt.savefig("sterile_numi.png")
  
  
  
  
  
