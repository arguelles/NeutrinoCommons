"""
Author  : C.A. Arguelles
Date    : 10/MAY/2011

Implements several bodies.

Log :
- Modified 11/FEB/2012 by C.A. Arguelles 
    + Created class to handle Earth with
    Cavity.
    + Create the innerclass class to handle
    inner class function calls.
- Modified 30/MAR/2012 by M.J. Bustamante
    + Created class to handle Stars genericaly
"""
import numpy as np
import scipy.interpolate as interpolate
import weakref, new
# my modules
import neutrinocommon
import neutrinocommon.tools.generaltools as gt
import neutrinocommon.physconst.physicsconstants as PC

# global variables
pc = PC.PhysicsConstants()

filepath = neutrinocommon.astro.__path__
datapath = filepath[0] + "/data/sun/"

class innerclass(object):
    """Descriptor for making inner classes.

    Adds a property 'owner' to the inner class, pointing to the outer
    owner instance.
    """

    # Use a weakref dict to memoise previous results so that
    # instance.Inner() always returns the same inner classobj.
    #
    def __init__(self, inner):
        self.inner= inner
        self.instances= weakref.WeakKeyDictionary()

    # Not thread-safe - consider adding a lock.
    #
    def __get__(self, instance, _):
        if instance is None:
            return self.inner
        if instance not in self.instances:
            self.instances[instance]= new.classobj(
                self.inner.__name__, (self.inner,), {'owner': instance}
            )
        return self.instances[instance]

class Vacuum():
    name = "Vacuum"
    Radius = 1.0

    # Vacuum region
    def density(self,track):
        return 0.0

    def rdensity(self,track):
        return 0.0

    def ye(self,track):
        return 1.0

    def RNC(self,track):
        return 1.0

    @innerclass
    class track():
        # Track information just for code compatibility
        def __init__(self,xini,xend):
            # x : current position [eV]
            # xini : initial position [eV]
            # xend : final position [eV]
            self.x = xini
            self.xini = xini
            self.xend = xend

class ConstantDensity():
    name = "Constant density"

    Radius = 1.0

    def __init__(self,rho,yye):
        self.rho = rho
        self.yye = yye

    # Vacuum region
    def density(self,track):
        return self.rho

    def ye(self,track):
        return self.yye

    def RNC(self,track):
        ye = self.yye
        return -0.5*(1.0-ye)/ye

    class track():
        # Track information just for code compatibility
        def __init__(self,xini,xend):
            # x : current position [eV]
            # xini : Initial position [eV]
            # xend : Final position [eV]
            self.x = xini
            self.xini = xini
            self.xend = xend

class VariableDensity():
    name = "Variable density"
    Radius = 1.0

    def __init__(self,rhomin,rhomax,yye):
        self.rhomin = rhomin
        self.rhomax = rhomax
        self.yye = yye

    def rdensity(self,r):
        # r = 1 => rhomax, r = 0 => rhomin
        m = (self.rhomax-self.rhomin)
        return self.rhomin + m*r

    def density(self,track):
        # Linear Density
        m = (self.rhomax-self.rhomin)/(track.xend-track.xini)
        return self.rhomin+m*(track.x-track.xini)

    def ye(self,track):
        # fix ye
        return self.yye

    def RNC(self,track):
        # fix ACC/ANC ratio
        ye = self.yye
        return -0.5*(1.0-ye)/ye
        
    @innerclass
    class track():
        def __init__(self,xini,xend):
            # x : current position [eV]
            # xini : Initial position [eV]
            # xend : Final position [eV]
            self.x = x
            self.xini = xini
            self.xend = xend

class Sun():
    name = "Sun"

    file = open(datapath+'bs05_agsop.dat','r')
    header, SSM = gt.hreadfilev2(file)

    ## Begin Loading Standar Solar Model ##
    SSM = np.array(SSM)
    M = SSM[:,0]
    R = SSM[:,1]
    T = SSM[:,2]
    Rho = SSM[:,3]
    P = SSM[:,4]
    Luminosity = SSM[:,5]
    X = SSM[:,6]
    ## End Loading Standard Solar Model ##

    SunRadius = 694439.0       # [km]
    Radius    = 694439.0       # [km]

    # Define interpolators
    # spline interpolator
    #inter_rdensity  = interpolate.UnivariateSpline(R,Rho)
    #inter_rxh       = interpolate.UnivariateSpline(R,X)

    # linear interpolator
    inter_rdensity  = interpolate.interp1d(R,Rho)
    inter_rxh       = interpolate.interp1d(R,X)

    def rdensity(self,r):
        # Density at a distance r from Sun center normalize to sun radius
        if(r<=self.R[-1] and r>= self.R[0]):
            return  self.inter_rdensity(r)
        elif(r<self.R[0]):
            return self.Rho[0]
        else :
            return 0.0

    def rxh(self,r):
        # mean molecular density at distance r from center
        if(r<=self.R[-1] and r>= self.R[0]):
            return  self.inter_rxh(r)
        elif(r<self.R[0]):
            return self.X[0]
        elif(r>self.R[-1]):
            return self.X[-1]
        else:
            return 0.0

    def rye(self,r):
        # mean molecular density at distance r from center
        return 0.5*(1.0+self.rxh(r))

    def rRNC(self,r):
        # ratio between neutral and charge current
        ye = self.rye(r)
        return -0.5*(1.0-ye)/ye

    def density(self,track):
        # Returns sun density at a given track position
        # considering radial motion
        # in this case geometry is trivial
        xkm = track.x/pc.km
        r = xkm/self.SunRadius
        if(r<=self.R[-1] and r>= self.R[0]):
            return np.max(self.inter_rdensity(r),0.0)
        elif(r<self.R[0]):
            return self.Rho[0]
        else :
            return 0.0
        #return self.rdensity(xkm/self.SunRadius)

    def ye(self,track):
        # Electron mean molecular weight
        # from J. Bacahall New Journal of Physics 6 (2004) 63, p.8
        xkm = track.x/pc.km
        r   = xkm/self.SunRadius
        return 0.5*(1.0+self.rxh(r))
        #return self.rye(xkm/self.SunRadius)

    def RNC(self,track):
        # ratio between neutral and charge current
        xkm = track.x/pc.km
        return self.rRNC(xkm/self.SunRadius)

    @innerclass
    class track():
        # in this case track information is trivial
        def __init__(self,xini,xend):
            # x : current position [eV]
            # xini : Initial radius [eV]
            # xendf : Final radius [eV]
            self.x = xini
            self.xini = xini
            self.xend = xend

class Star():
    """
    Defining a Star
    """
    name = "Star"
    
    def __init__(self,Radius,InitialRadius,InitialDensity):
        self.Radius = Radius
        self.InitialRadius = InitialRadius
        self.InitialDensity = InitialDensity

    def rdensity(self,r):
        # Electron Density at a distance r from Sun center normalize to sun radius
        if(r<=self.Radius and r>= self.InitialRadius):
            return self.InitialDensity*np.exp( - r/self.InitialRadius)
        elif(r<self.InitialRadius):
            return self.InitialDensity
        else :
            return 0.0

    def density(self,track):
        # Returns star density at a given track position
        # considering radial motion
        # in this case geometry is trivial
        xkm = track.x/pc.km
        return self.rdensity(xkm)
        
    def rye(self,r):
        # electron mean molecular density at distance r from center
        return 0.86#0.5*(1.0+self.rxh(r))

    def ye(self,track):
        # Electron mean molecular weight
        # from J. Bacahall New Journal of Physics 6 (2004) 63, p.8
        xkm = track.x/pc.km
        return self.rye(xkm)

    @innerclass
    class track():
        # in this case track information is trivial
        def __init__(self,xini,xend):
            # x : current position [eV]
            # xini : Initial radius [eV]
            # xendf : Final radius [eV]
            self.x = xini
            self.xini = xini
            self.xend = xend

class Earth():
    name = "Earth"

    EarthRadius = 6371.0    #[km]
    Radius = 6371.0    #[km]

    def rdensity(self,x):
        # Calcula la densidad de la Tierra segun el PREM
        # R. Gandhi et al. Astroparticle Phys. 5, 81-110 (1996)
        # Arxiv : 9512364 pag. 23
        # x is adimentional radius : x = 0 : center, x = 1 : Earthradius
        r = self.EarthRadius*x
        if r <= 1221.50 :
            dne = 13.08850-8.83810*x**2
        elif r>=1221.50 and r<3480 :
            dne=12.58150-1.26380*x-3.64260*x**2.-5.5280*x**3.
        elif r >=3480.0 and r < 5701.0 :
            dne=7.95650-6.47610*x+5.52830*x**2.-3.08070*x**3.
        elif r >= 5701.0 and r<5771.0 :
            dne=5.31970-1.48360*x
        elif r>=5771.0 and r<5971.0 :
            dne=11.24940-8.02980*x
        elif r>=5971.0 and r<6151.0 :
            dne=7.10890-3.80450*x
        elif r>=6151.0 and r<6346.60 :
            dne=2.6910+0.69240*x
        elif r >= 6346.60 and r < 6356.0 :
            dne = 2.9
        elif r >= 6356.0 and r < 6368 :
            dne = 2.6
        elif r<= self.EarthRadius :
            dne = 1.020
        elif r>=self.EarthRadius :
            dne=0.0
        return dne

    def density(self,track):
        # Returns earth density at a given track position
        xkm = track.x/pc.km
        # patch to conv. L back to km
        r = np.sqrt(self.EarthRadius**2-(track.L/pc.km)*xkm+xkm*xkm)
        return self.rdensity(r/self.EarthRadius)

    def ye(self,track):
        # Electron mean molecular weight
        return 0.494

    @innerclass
    class track():
        x = float
        def __init__(self,xini,xend,L):
            # L : baseline [eV^-1]
            # x : current position [eV^-1]
            # xini : initial position [eV^-1]
            # xend : final position [eV^-1]
            self.L = L
            self.x = xini
            self.xini = xini
            self.xend = xend

class EarthWithCavity():
    name = "Earth with cavity"

    Radius = 6371.0*pc.km
    EarthRadius = Radius/pc.km

    def __init__(self,rho_cavity,ye_cavity,wx,wy,th_rot,deep_cavity,th_cavity):
        """ Initializes EarthWithCavity object by defining the cavities
        properties and position.

        @type   rho_cavity  :   float
        @param  rho_cavity  :   cavity density [gr/cm^3]
        @type   ye_cavity   :   float
        @param  ye_cavity   :   electron fraction [dimensionless]
        @type   wx          :   float
        @param  wx          :   cavity horizontal axis (x-width) [eV^-1]
        @type   wy          :   float
        @param  wy          :   cavity vertical axis (y-width) [eV^-1]
        @type   th_rot      :   float
        @param  th_rot      :   axis orientations with respect to the equator [rad]
        @type   deep_cavity :   float
        @param  deep_cavity :   distance of cavity center to Earth's surface measure along the vertical axis (y of cavity center) [eV^-1]
        @type   th_cavity   :   float
        @param  th_cavity   :   angle measure from the meridian (+y) to the cavity center [rad]
        """
        # cavity matter properties
        self.rho_cavity = rho_cavity
        self.ye_cavity = ye_cavity
        # cavity geometry
        self.wx = wx
        self.wy = wy
        self.th_rot = th_rot
        self.deep_cavity = deep_cavity
        self.th_cavity = th_cavity

    def rdensity(self,x):
        # Calcula la densidad de la Tierra segun el PREM
        # R. Gandhi et al. Astroparticle Phys. 5, 81-110 (1996)
        # Arxiv : 9512364 pag. 23
        # x is adimentional radius : x = 0 : center, x = 1 : Earthradius
        EarthRadius = self.EarthRadius

        r = EarthRadius*x

        if r>= EarthRadius:
            return 0.0
        elif r>= 6369.0 and r< EarthRadius:
            return 1.020
        elif r >= 6356.0 and r < 6368.0 :
            return 2.6
        elif r >= 6346.60 and r < 6356.0 :
            return 2.9
        elif r>=6151.0 and r<6346.60 :
            return 2.6910+0.69240*x
        elif r>=5971.0 and r<6151.0 :
            return 7.10890-3.80450*x
        elif r>=5771.0 and r<5971.0 :
            return 11.24940-8.02980*x
        elif r >= 5701.0 and r<5771.0 :
            return 5.31970-1.48360*x
        elif r >=3480.0 and r < 5701.0 :
            return 7.95650-6.47610*x+5.52830*x**2.-3.08070*x**3.
        elif r>=1221.50 and r<3480 :
            return 12.58150-1.26380*x-3.64260*x**2.-5.5280*x**3.
        elif r <= 1221.50 :
            return 13.08850-8.83810*x**2
        else:
            return 0.0

    def density(self,track):
        """ Return the density at a given track position.

        @type   track  :   track
        @param  track  :   track object

        @rtype         :   float
        @return        :   density [gr/cm^3]
        """
#        if self.isInCavity(track) :
        if track.x >= track.xin and track.x <= track.xout :
            return self.rho_cavity
        else :
            xkm = track.x/pc.km
            # patch to conv. L back to km
            r = np.sqrt(self.EarthRadius**2-(track.L/pc.km)*xkm+xkm*xkm)
            return self.rdensity(r/self.EarthRadius)

    def ye(self,track):
        """ Return the electron fraction at a given track position.

        @type   track  :   track
        @param  track  :   track object

        @rtype         :   float
        @return        :   electron fraction [dimensionless]
        """

#        if self.isInCavity(track) :
        if track.x >= track.xin and track.x <= track.xout :
            return self.ye_cavity
        else :
            # assuming crust value
            return 0.494

    def isInCavity(self,track):
        """ Return True if the particle is inside of the cavity.

        @type   track  :   track
        @param  track  :   track object

        @rtype              :   boolean
        @return             :   True : inside cavity, False : outside cavity
        """
        # We'll do this calculation in a coordinate system
        # centered in the Earth whose x-axis is align with
        # the Earth's equator. In this system (0,0) is the
        # Earth's center.
        Earth_radius = self.Radius

        # neutrino source xy-coordinates
        # minus sign due to the definition of track.phi
        x_source = -Earth_radius*np.sin(track.phi)
        y_source = Earth_radius*np.cos(track.phi)

        # neutrino xy-coordinates
        xcurrent = track.x
        x_neutrino = x_source + xcurrent*np.cos(track.th)
        y_neutrino = y_source + xcurrent*np.sin(track.th)

        # cavity center xy-coordinates
        R_cavity = (Earth_radius - self.deep_cavity)
        # minus sign due to the definition of th_cavity
        #x_cavity = -R_cavity*np.sin(self.th_cavity)
        #y_cavity = R_cavity*np.cos(self.th_cavity)
        # with a "horizontal" definition of deep
        x_cavity = -R_cavity*np.tan(self.th_cavity)
        y_cavity = R_cavity

        # we need to check if the (x_neutrino,y_neutrino)
        # is in the cavity ellipse

        # considering the cavity rotation
        th_rot = self.th_rot
        x_rot = x_neutrino*np.cos(th_rot) - y_neutrino*np.sin(th_rot)
        y_rot = x_neutrino*np.sin(th_rot) + y_neutrino*np.cos(th_rot)

        # a and b of the ellipse
        a_cav = self.wx/2.0
        b_cav = self.wy/2.0

        return ((x_rot-x_cavity)/a_cav)**2+((y_rot-y_cavity)/b_cav)**2 <= 1.0

    def GetBaseline(self,track):
        """ Return neutrino baseline in kilometers.

        @type   self  :   track
        @param  self  :   track object

        @rtype        :   float
        @return       :   baseline [km]
        """
        # the angle that oposes the baseline is
        baseline_angle = 2.0*(track.phi-track.th)
        # using the law of cosines
        return np.sqrt((2.0*self.Radius**2)*(1.0-np.cos(baseline_angle)))

    def GetTransversedWidthAndDistance(self,track):
        """ Return neutrino traverse cavity width
        in kilometers. Also returns the distance
        to the cavity border in kilometers.

        @type   self  :   track
        @param  self  :   track object

        @rtype        :   float
        @return       :   (traverse cavity width [km], distance to border [km])
        """
        if self.wx != 0.0 and self.wy != 0.0 :
            Earth_radius = self.Radius

            # naming angles
            phi = track.phi
            th = track.th
            th_cavity = self.th_cavity
            rth = self.th_rot

            # neutrino source xy-coordinates
            # minus sign due to the definition of track.phi
            xs = -Earth_radius*np.sin(phi)
            ys = Earth_radius*np.cos(phi)

            # cavity center xy-coordinates
            R_cavity = (Earth_radius - self.deep_cavity)
            # minus sign due to the definition of th_cavity
            # with a "radial" definition of the deep
            #xc = -R_cavity*np.sin(th_cavity)
            #yc = R_cavity*np.cos(th_cavity)
            # with a "horizontal" definition of deep
            xc = -R_cavity*np.tan(th_cavity)
            yc = R_cavity

            # cavity size renaming
            wx = self.wx/2.0
            wy = self.wy/2.0

            # slope
            m = np.tan(th)

            # AUXILIARY COEFFICIENTS
            # rotated ellipse : X^2cxx+Y^2yxx++Xcx+Ycy+CC = 0
            cxx = np.cos(rth)**2/wx**2 + np.sin(rth)**2/wy**2
            cyy = np.sin(rth)**2/wx**2 + np.cos(rth)**2/wy**2
            cxy = np.sin(2*rth)*(1.0/(wy*wy)-1.0/(wx*wx))
            #cx = -2*(xc*np.cos(rth)/(wx*wx)+yc*np.sin(rth)/(wy*wy))
            #cy = -2*(yc*np.cos(rth)/(wy*wy)-xc*np.sin(rth)/(wx*wx))
            #cc = (xc/wx)**2+(yc/wy)**2-1.0
            # equivalent to previous lines (OK - TESTED)
            cx = -(2*cxx*xc+cxy*yc)
            cy = -(2*cyy*yc+cxy*xc)
            cc = cxx*xc**2+cyy*yc**2+xc*yc*cxy-1.0
            # CUADRATIC COEFFICIENTS ASUMING Ax^2 + Bx + C = 0
            a2 = m*m*cyy+m*cxy+cxx
            b2 = 2*m*(ys-xs*m)*cyy+(ys-xs*m)*cxy+m*cy+cx
            c2 = (ys-xs*m)**2*cyy+cy*(ys-xs*m)+cc
            delta2 = b2*b2-4*a2*c2
            if delta2 <= 0.0 :
                # beam doesn't pass through the cavity
                return -1,-1
            else :
                # CAVITY INTERSECTION COORDINATES
                c1x = -b2/(2.0*a2)+np.sqrt(delta2)/(2.0*a2)
                c1y = ys + m*(c1x-xs)
                c2x = -b2/(2.0*a2)-np.sqrt(delta2)/(2.0*a2)
                c2y = ys + m*(c2x-xs)
                # CALCULATING NEW WIDTH
                w2 = np.sqrt((c1x-c2x)**2+(c1y-c2y)**2)
                # CALCULATING NEW DISTANCE TO CAVITY
                d2 = np.sqrt((c2x-xs)**2+(c2y-ys)**2)
                return w2/pc.km,d2/pc.km
        else :
            return -1,-1
            
    @innerclass
    class track():
        def __init__(self,phi,th):
            """ Creates a track object related
            to this body. Also calculates the
            baseline and stores it.

            @type   phi  :   float
            @param  phi  :   specifies the source position measure from (+y) counterclockwise [rad].
            @type   th   :   float
            @param  th   :   angle between the baseline and the equator line [rad].

            @rtype        :   track
            @return       :   track object
            """
            ### USUAL ###
            # x : current position [eV^-1]
            # xini : initial position [eV^-1]
            # xend : final position [eV^-1]
            # L : baseline [eV^-1]
            Radius = 6371.0*pc.km

            # set baseline [eV^-1]
            self.phi = phi
            self.th = th

            # the angle that oposes the baseline is
            baseline_angle = 2.0*(phi-th)
            # using the law of cosines
            self.L = np.sqrt((2.0*Radius**2)*(1.0-np.cos(baseline_angle)))

            self.x = 0.0
            self.xini = 0.0
            self.xend = self.L

            (w,d) = EarthWithCavity.GetTransversedWidthAndDistance(self.owner,self)
            self.xin = d*pc.km
            self.xout = (w+d)*pc.km

class EarthZenith():
    name = "Earth"

    Radius = 6371.0*pc.km
    EarthRadius = Radius/pc.km

    def __init__(self):
        """ Initializes EarthWithCavity object by defining the cavities
        properties and position.
        """
        # nothing to initialize

    def rdensity(self,x):
        # Calcula la densidad de la Tierra segun el PREM
        # R. Gandhi et al. Astroparticle Phys. 5, 81-110 (1996)
        # Arxiv : 9512364 pag. 23
        # x is adimentional radius : x = 0 : center, x = 1 : Earthradius
        EarthRadius = self.EarthRadius

        r = EarthRadius*x

        if r>= EarthRadius:
            return 0.0
        elif r>= 6369.0 and r< EarthRadius:
            return 1.020
        elif r >= 6356.0 and r < 6368.0 :
            return 2.6
        elif r >= 6346.60 and r < 6356.0 :
            return 2.9
        elif r>=6151.0 and r<6346.60 :
            return 2.6910+0.69240*x
        elif r>=5971.0 and r<6151.0 :
            return 7.10890-3.80450*x
        elif r>=5771.0 and r<5971.0 :
            return 11.24940-8.02980*x
        elif r >= 5701.0 and r<5771.0 :
            return 5.31970-1.48360*x
        elif r >=3480.0 and r < 5701.0 :
            return 7.95650-6.47610*x+5.52830*x**2.-3.08070*x**3.
        elif r>=1221.50 and r<3480 :
            return 12.58150-1.26380*x-3.64260*x**2.-5.5280*x**3.
        elif r <= 1221.50 :
            return 13.08850-8.83810*x**2
        else:
            return 0.0

    def density(self,track):
        """ Return the density at a given track position.

        @type   track  :   track
        @param  track  :   track object

        @rtype         :   float
        @return        :   density [gr/cm^3]
        """
        xkm = track.x/pc.km
        # patch to conv. L back to km
        r = np.sqrt(self.EarthRadius**2-(track.L/pc.km)*xkm+xkm*xkm)
        return self.rdensity(r/self.EarthRadius)

    def ye(self,track):
        """ Return the electron fraction at a given track position.

        @type   track  :   track
        @param  track  :   track object

        @rtype         :   float
        @return        :   electron fraction [dimensionless]
        """
        # need to change
        return 0.494
        
    @innerclass
    class track():
        def __init__(self,phi):
            """ Creates a track object related
            to this body. Also calculates the
            baseline and stores it.

            @type   phi  :   float
            @param  phi  :   zenith angle [rad].

            @rtype        :   track
            @return       :   track object
            """
            ### USUAL ###
            # x : current position [eV^-1]
            # xini : initial position [eV^-1]
            # xend : final position [eV^-1]
            # L : baseline [eV^-1]
            Radius = 6371.0*pc.km

            # set baseline [eV^-1]
            self.phi = phi
            self.th = th

            # the angle that oposes the baseline is
            baseline_angle = 2.0*(phi-th)
            # using the law of cosines
            self.L = np.sqrt((2.0*Radius**2)*(1.0-np.cos(baseline_angle)))

            self.x = 0.0
            self.xini = 0.0
            self.xend = self.L

            (w,d) = EarthWithCavity.GetTransversedWidthAndDistance(self.owner,self)
            self.xin = d*pc.km
            self.xout = (w+d)*pc.km

class AGN():
    name = "AGN"

    Radius = 0.02*3.08568e13    #[km]

    def __init__(self,rhomax,yye):
        self.rhomax = rhomax # [gr/cm^3]
        self.yye = yye       # [dimensionless]

    def rdensity(self,x):
        # Estimates the AGN density using an exponential profile.
        # x is adimentional radius : x = 0 : center, x = 1 : AGN radius
        return self.rhomax*np.exp(-x*10.0)

    def density(self,track):
        # Returns earth density at a given track position
        xparsec = track.x/pc.km
        return self.rdensity(xparsec/self.Radius)

    def ye(self,track):
        # Electron mean molecular weight
        return self.yye

    @innerclass
    class track():
        x = float
        def __init__(self,xini,xend):
            # x : current position [eV]
            # xini : initial position [eV]
            # xend : final position [eV]
            self.x = xini
            self.xini = xini
            self.xend = xend
