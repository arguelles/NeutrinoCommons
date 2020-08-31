===============
NeutrinoCommons
===============

NeutrinoCommons provides a various python scripts to perform neutrino
oscillation experiments. Here in you will find a neutrino propagation
algorithm as well as neutrino cross section data, several neutrino 
sources (e.g. Sun, Earth, Accelerators, etc.).

Critical subroutines are implemented in cython, therefore providing
considerable speed gain. It is most useful for task involving the 
calculation of neutrino oscillations and also estimating event numbers.
Typical usage looks like this::

    #!/usr/bin/env python

    import neutrinocommon.astro.body as body
    import neutrinocommon.neu.neuosc as neuosc
    import neutrinocommon.physconst.physicsconstants as PC
    
    pc = PC.PhysicsConstants()
    
    # Defines the body where the neutrino will propagate.
    earth = body.Earth()
    # Defines the trajectory within the body that the neutrino will take.
    baseline = 1000.0*pc.km
    xend = baseline
    xini = 0.0

    track = earth.track(xini,xend,baseline)    
    
    The oscillation probability is implemented by the two functions
    CalNeuOscPy and CalNeuOscGSL. The later is much more faster and complete
    code, but requieres the GSL-python libraries, while the first one only requieres
    scipy and numpy.

    Since this functions are usually called several times they requiere the 
    mass square difference in the flavor basis

    We create this matrix for this parameter set
    flavorM2 = no.flavorM2(pc)

    We also need to specify the initial neutrino flavor and its energy, and also
    Set if we are talking about a neutrino or an antineutrino. 

    E = 80.0*pc.MeV
    ineu = 0
    pc.neutype = "antineutrino"

    print "Oscillation probabilities : "
    print no.CalNeuOscPy(ineu,E,track,body,flavorM2,pc)
    print no.CalNeuOscGSL(ineu,E,track,body,flavorM2,pc)
    
    
Please view the demo.py in the DemoProject.

==================
Instalations Notes
==================

- Python Requierments:
    + numpy
    + scipy
    + datetime
    + os
    + fnmatch
    + pyGSL
    + matplotlib

- PyGSL Requiers :
    + pythonX.Y-dev
    + gls-config
    + libgsl0-dev
    
- Other Requierments:
    + cython
    + pyrex

