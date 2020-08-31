"""
Generate wrapper for TAUOLA - simple tau decay. 

Run : python setup.py build_ext --inplace
"""
# local variables
# this variables should be modified accordingly
# MINIMUN  :  HEPMC, TAUOLA
# OPTIONAL : MCTESTER, ROOT, PYTHIA8

# home-root
# in linux-like systems
#home = "/home/carguelles/Programs/"
# in MACOSX
home = "/Users/carguelles/Programs/"
rel_dir = "../../../linked/"

# specify location of pythia
#pythialib = home+"pythia8/pythia8170/lib/archive"
#pythiainc = home+"pythia8/pythia8170/include"
# libraries : hepmcinterface

# specify location of ROOT
#rootlib = home+"i3/ports/root-v5.27.00/lib"
#rootinc = home+"i3/ports/root-v5.27.00/include"
#rootlib = home+"root/build/lib"
#rootinc = home+"root/build/include"


# specify location of MCTESTER
#MCTESTERlib = home+"mctester/MC-TESTER/lib"
#MCTESTERinc = home+"mctester/MC-TESTER/include"
# libraries : HepMCEvent,HEPEvent,MC-Tester

# specify location of HepMC
#HepMClib = home+"hepmc/build/lib"
HepMClib = rel_dir+"hepmc/build/lib"

# specify location of TAUOLA
#TAUOLAlib = home+"tauola/TAUOLA/lib"
#TAUOLAinc = home+"tauola/TAUOLA/include/Tauola"
TAUOLAlib = rel_dir+"tauola/src/lib"
TAUOLAinc = rel_dir+"tauola/src/include/Tauola"

# specify local libraries
locallib = "/usr/local/lib"
localinc = "/usr/local/include"

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ROOTLIBS = "Gui","Gpad","Hist","Graf","Graf3d","Tree","Rint","Postscript","Matrix","Physics","MathCore","RIO","Net","Thread","Core","Cint","m","dl"
PYTHIALIBS = "pythia8","hepmcinterface"

ext_modules = [Extension("taudecay", ["taudecay.pyx"],language="c++",include_dirs=[localinc,TAUOLAinc],libraries=["HepMC","TauolaCxxInterface"],extra_link_args=["-fPIC","-rdynamic","-L"+TAUOLAlib,"-L"+HepMClib])]

#ext_modules = [Extension("taudecay", ["taudecay.pyx"],language="c++",include_dirs=[localinc,pythiainc,rootinc,TAUOLAinc],libraries=["HepMC","TauolaCxxInterface","lhapdfdummy"],extra_link_args=["-fPIC","-rdynamic","-L"+rootlib,"-L"+pythialib,"-L"+TAUOLAlib,"-L"+HepMClib])]

setup(
  name = 'taudecay',
  cmdclass = {'build_ext': build_ext},
  ext_modules = ext_modules
)
