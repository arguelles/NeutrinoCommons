"""
Author  : C.A. Arguelles
Date    : 10/MAY/2011

This script constains the MINOS experiment information
and configuration. It is organized in several classes
for its easier setup.

Log:
- Modified on 23/FEB/2012 by C.Arguelles
    + Review the code to make it compatible with the neutrino
    commons library.
"""

## NOTE ##
# This enthought trait libraries should be reconsidered
# to use only python standard functions. Or consider the 
# option of including this in the setup file.

from enthought.traits.api import *
from enthought.traits.ui.api import *
from matplotlib.figure import Figure
from threading import Thread
# my modules
import neutrinocommon.physconst.physicsconstants as PC

class MINOSConfiguration(HasTraits):
  """
  Minos plotter configuration
  """
  
  # PLOT PARAMETERS
  th_23_step = CFloat(0.01,
              desc = "th_23_step plotting step.",
              label = "Step sin2",)
  th_23_min = CFloat(0.45,#0.0,
              desc = "minimum value of th_23.",
              label = "Minimum sin2",)
  th_23_max = CFloat(0.79,
              desc = "maximum value of th_23.",
              label = "Maximum sin2",)
  dm32sq_step = CFloat(0.5e-4,
              desc = "dm^2 plotting step.",
              label = "Step dm2",)
  dm32sq_min = CFloat(1.5e-3,
              desc = "minimum value of dm^2.",
              label = "Minimum dm2",)
  dm32sq_max = CFloat(3.5e-3,
              desc = "maximum value of dm^2.",
              label = "Maximum dm2",)
  
  #ranged = CBool(True,
              #label = "Ranged")
  #plot = CBool(True,
              #label = "Plot")
  notemp = CBool(False,
              label = "Store temporal file")
  optimize = CBool(True,
              label = "Use optimization with external file")
  usedata = CBool(False,
              label = "Use external data to make plots")
  experiment = CBool(False,
              label = "Generate experiment beam specter")
  format = CStr("png",
              label = "Image format")

class MINOSExperiment(MINOSConfiguration,PC.PhysicsConstants):
  """
  Minos Experiment Object
  """
  # BASELINE DISTANCE
  baseline = CFloat(735.0,
              desc = "Distance from Far Detector [km]",
              label= "Baseline length",)
  
  # ENERGY PARAMETERS
  minenergy = CFloat(0.01,
              desc = "Minimum energy in the integration range [GeV].",
              label = "Minimum energy",)
  maxenergy = CFloat(50.0,
              desc = "Maximum energy in the integration range [GeV].",
              label = "Maximum energy",)
  binsize = CFloat(1.0,
              desc = "Energy binsize [GeV].",
              label = "Energy binsize",)

  # NEUTRINO TYPE
  neutype = CFloat(0,
              desc = "Neutrino Type. 0 : neutrino, 1 : antineutrino.",
              label = "Neutrino type",)
  
  # BEAM
  class Beam():
    POT_neutrino 	= 7.20e20	# Neutrino run Protons-on-target
    POT_antineutrino 	= 1.71e20	# Antineutrino run Protons-on-target
  
  # FAR DETECTOR
  class FarDet():
    # From MINOS document 6565-v2 p. 165
    # view : http://arxiv.org/pdf/0805.3170v2
    mass = 5.4e6#2673.0*1.0e3#3.3*1.0e6	#	[kg]
    distance = 735.34			#	[km]
    
    class Material():
      name = 'Fe'
      NucleonNum = 56.0			#	[nucleons/nuclei]
      AtomicWeight = 55.845		#	[gr/mol]
      
  class NearDet():
    # From MINOS document 6565-v2 p. 165
    distance = 0.84			#	[kg]
  
  # OSCILATION CONFIGURATION
  class NeuOscConf():
    generations	= CBool(True,
		desc = "If True : 3 generation formalism will be used. Else two generation formalism will be used.",
		label = "Neutrino oscillations")
    oscillation = CBool(True,
		desc = "If True : neutrino oscillations will be considered.",
		label = "Neutrino oscillations")
    matter 	= CBool(True,
		desc = "If True : matter effects will be taken into account. Else vacuum oscillations will be used.",
		label = "Matter")
              
  # REAL OSCILLATION VALUES FOR NEUTRINO
  if neutype == 0 :
    sin2real = CFloat(1.0,
	 desc = "Real sin(2\theta)^2.",
	 label= "Real sin2",)
    dm2real = CFloat(2.35e-3,
	 desc = "Real dm23^2.",
	 label= "Real dm2",)
  elif neutype == 1:
    sin2real = CFloat(0.86,
	 desc = "Real sin(2\theta)^2.",
	 label= "Real sin2",)
    dm2real = CFloat(3.36e-3,
	 desc = "Real dm23^2.",
	 label= "Real dm2",)
  # OTHER
  globalfact = CFloat(1.0,
	desc = "arbitrary event-number multiplicative factor",
	label = "Global multiplative factor")
