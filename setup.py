from distutils.core import setup
from distutils.extension import Extension

import numpy as np

## compiling cython extensions

try:
    from Cython.Distutils import build_ext
except ImportError:
    use_cython = False
else:
    use_cython = True
    
cmdclass = {}
ext_modules = []

if use_cython :
    ext_modules += [
        Extension('neutrinocommon.shareobjects.optgeneraltools',["neutrinocommon/shareobjects/cython_codes/optgeneraltools/optgeneraltools.pyx"],extra_compile_args=['-I'+np.__path__[0]+'/core/include']),        
        Extension('neutrinocommon.shareobjects.optneuosc',["neutrinocommon/shareobjects/cython_codes/optneuosc/optneuosc.pyx"],extra_compile_args=['-I'+np.__path__[0]+'/core/include']),

    ]
    cmdclass.update({ 'build_ext': build_ext })
else:
    ext_modules += [
        Extension('neutrinocommon.shareobjects.optgeneraltools',["neutrinocommon/shareobjects/cython_codes/optgeneraltools/optgeneraltools.c"],extra_compile_args=['-I'+np.__path__[0]+'/core/include']),        
        Extension("neutrinocommon.shareobjects.optneuosc", ["neutrinocommon/shareobjects/cython_codes/optneuosc/optneuosc.c"],extra_compile_args=['-I'+np.__path__[0]+'/core/include']),
    ]

## setting up neutrinocommons

setup(
    name='neutrinocommon',
    version='0.3.0',
    author='C.A.Arguelles',
    author_email='carguellesdel@gmail.com',
    packages=['neutrinocommon','neutrinocommon.astro', 'neutrinocommon.exp','neutrinocommon.exp.icecube','neutrinocommon.exp.minos','neutrinocommon.exp.miniboone','neutrinocommon.neu','neutrinocommon.physconst','neutrinocommon.plotters','neutrinocommon.shareobjects','neutrinocommon.tools','neutrinocommon.beams'],
    package_data={'neutrinocommon.astro': ['data/dm/DMnuProdFluxes/*.dat','data/dm/DMnuProdParameters/*.dat','data/dm/wimpsim/*.dat','data/sun/*.dat'],'neutrinocommmon.exp':['icecube/data/*.dat','miniboone/data/*.dat','minos/data/*.dat'],'neutrinocommon.neu':['data/xs/*.dat','bin/xs/*.exe','bin/osc/*.exe'],'neutrinocommon.shareobjects':['*.so']},
    description='Several tools for neutrino phenomenology elaborated by the HEP-PUCP theory group.',
    long_description=open('README.txt').read(),
    cmdclass = cmdclass,
    ext_modules = ext_modules
)

## compiling fortran extensions

from numpy.distutils.core import setup, Extension

setup(
ext_modules = [Extension('neutrinocommon.shareobjects.nudsde',["neutrinocommon/shareobjects/cython_codes/nudsde/nudsde_python.f"],include_dirs=['neutrinocommon/linked/nusigma/inc'],library_dirs=['neutrinocommon/linked/nusigma/lib'],libraries=['nusigma','xcm','pythia'])],
)

