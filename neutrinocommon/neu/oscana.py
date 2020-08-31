""" 
Author  : C.A. Arguelles
Date    : 10/MAY/2011

This package contains formulaes to evaluate the neutrino
oscillation probability analytical, in cases such as, constant
density, vacuum, or others. Also the decoherence formalism
is implemented.

Log :
- Modified on 23/ABR/2012 by C.A.Arguelles
    + Adapted the code to work with the neutrino
    commons package.
- Modified on 31/MAY/2012 by C.A.Arguelles
    + Added the adiabatic formulas.
"""

# standard python modules
import numpy as np
import subprocess
from scipy import integrate
# my modules
import neutrinocommon
import neutrinocommon.neu.neuosc as no
import neutrinocommon.tools.generaltools as gt
import neutrinocommon.physconst.physicsconstants as PC

# global variables
pc = PC.PhysicsConstants()

#===============================================================================
# BEGIN DECOHERENCE FORMULAE
#===============================================================================

"""
References for this section:

1) A Study on quantum decoherence phenomena with three generations of neutrinos
A.M. Gago, E.M. Santos, W.J.C. Teves, R. Zukanovich Funchal
arXiv : hep-ph/0208166

2) Quantum decoherence and neutrino data
G. Barenboim,N.E. Mavromatos, A. Waldron-Lauda
arXiv : hep-ph/0603028
"""

def RhoDecoherence(alpha,param):
	"""
    Ref. (1). Appendix B : coeficients from Apendix B, ec. B4.
	
    @type   alpha     		:      integer
    @param  alpha    		:      matrix index
    @type   param     		:      physics constants
    @param  param   		:      physics constants
    
	"""
	U = no.mixmatrix(param).U
	UCT =  no.mixmatrix(param).UCT
	
	Rho = np.zeros([9],complex)
	
	Rho[0] = np.sqrt(2.0/3.0)
	Rho[1] = 2.0*np.real(UCT[0,alpha]*U[alpha,1])
	Rho[2] = -2.0*np.imag(UCT[0,alpha]*U[alpha,1])
	Rho[3] = gt.norm(U[alpha,0])**2 - gt.norm(U[alpha,1])**2
	Rho[4] = 2.0*np.real(UCT[0,alpha]*U[alpha,2])
	Rho[5] = -2.0*np.imag(UCT[0,alpha]*U[alpha,2])
	Rho[6] = 2.0*np.real(UCT[1,alpha]*U[alpha,2])
	Rho[7] = -2.0*np.imag(UCT[1,alpha]*U[alpha,2])
	Rho[8] = (1.0/np.sqrt(3.0))*(gt.norm(U[alpha,0])**2+gt.norm(U[alpha,1])**2-2.0*gt.norm(U[alpha,2])**2)
	
	return Rho

def Omega(E,param):
	"""
	# A Study on quantum decoherence phenomena with three generations of neutrinos
	# Coeficients from Apendix B, ec. B4
	"""
	# assuming p aprox= E 
	omega = np.zeros(4,complex)
	omega[0] = np.sqrt(complex((param.gamma[1]-param.gamma[0])**2-param.dm21sq**2/E**2,0.0))						# In literature Omega12
	omega[1] = np.sqrt(complex((param.gamma[4]-param.gamma[3])**2-param.dm31sq**2/E**2,0.0))						# In literature Omega13
	omega[2] = np.sqrt(complex((param.gamma[6]-param.gamma[5])**2-param.dm32sq**2/E**2,0.0))						# In literature Omega23
	omega[3] = np.sqrt(complex((param.gamma[2]-param.gamma[7])**2,0.0))												# In literature Omega38 when D38 = 0.0
	#omega[3] = np.sqrt(complex((param.gamma[2]-param.gamma[7])**2+4.0*(param.gamma[2]-param.gamma[7])**2,0.0))		# In literature Omega38 -> wrong!!!
	return omega

def ProbDecoherence(a,b,E,L,param):
	"""
	# Arxiv : 0603028 - Barenboim // from alberto's paper
	# Probability of neutrino oscillation + decoherence in 3g
	# assuming t aprox= L
	# P(neu_a -> neu_b)
	#	E : Energy 	 [eV]
	#	L : Distance [eV^-1]
	# Setting Omega38 = 0.0
	"""
	param.Refresh()
	ra = RhoDecoherence(a, param)
	rb = RhoDecoherence(b, param)
	omg = Omega(E,param)
	gam = param.gamma

	pab = (1.0/3.0)
	pab = pab + 0.5*((ra[1]*rb[1]+ra[2]*rb[2])*np.cos(gt.norm(omg[0])*L/2.0)+((2.0*param.dm21sq*(ra[1]*rb[2]-ra[2]*rb[2])+(gam[1]-gam[0])*(ra[1]*rb[1]-ra[2]*rb[2]))/gt.norm(omg[0]))*np.sin(gt.norm(omg[0])*L/2.0))*np.exp(-L*(gam[0]+gam[1])/2.0)
	pab = pab + 0.5*((ra[4]*rb[4]+ra[5]*rb[5])*np.cos(gt.norm(omg[1])*L/2.0)+((2.0*param.dm31sq*(ra[4]*rb[5]-ra[5]*rb[4])+(gam[3]-gam[4])*(ra[3]*rb[3]-ra[4]*rb[4]))/gt.norm(omg[1]))*np.sin(gt.norm(omg[1])*L/2.0))*np.exp(-L*(gam[3]+gam[4])/2.0)
	pab = pab + 0.5*((ra[6]*rb[6]+ra[7]*rb[7])*np.cos(gt.norm(omg[2])*L/2.0)+((2.0*param.dm32sq*(ra[6]*rb[7]-ra[7]*rb[6])+(gam[5]-gam[6])*(ra[5]*rb[5]-ra[6]*rb[6]))/gt.norm(omg[2]))*np.sin(gt.norm(omg[2])*L/2.0))*np.exp(-L*(gam[5]+gam[6])/2.0)
	pab = pab + 0.5*( ra[3]*rb[3]*np.exp(-gam[2]*L) + ra[8]*rb[8]*np.exp(-gam[7]*L))
	
	return pab.real

def PmuSurvivalProbDecoherence(E,L,param):
	"""
	# Arxiv : Search for quantum gravity with IceCube and high energy atmospheric neutrinos
	# ec (1)
	"""
	gam = param.gamma
	m = np.abs((param.gamma[6]-param.gamma[5])**2-param.dm32sq**2/E**2)
	pmumu = (1.0/3.0)
	pmumu = pmumu + 0.5*(np.exp(-gam[2]*L)*np.cos(param.th23)**4+(1.0/12.0)*np.exp(-gam[7]*L)*(1.0-3.0*np.cos(2.0*param.th23))**2)
	pmumu = pmumu + 0.5*(4.0*np.exp(-L*(gam[5]+gam[6])/2.0))*np.cos(param.th23)**2*np.sin(param.th23)**2*(np.cos(L*np.sqrt(m)/2.0)+np.sin(L*np.sqrt(m)/2.0)*(gam[5]-gam[6])/np.sqrt(m))
	return pmumu

def AvgDecoherenceProbability(a,b,alpha,param):
	""" Calculate the average decoherence probabily with a given set of std osc. parameters
	considering the oscillation exponent * distance as a single variable.

        @type   a     		:      integer
        @param  a    		:      initial neutrino flavor
	@type   b     		:      integer
        @param  b     		:      final neutrino flavor
	@type   alpha	    	:      float
        @param  alpha     	:      decoherence exponent = gamma*L
        @type   param           :      physicsconstant
        @param  param           :      set of physical parameters use for estimating the probability.
        
        @rtype                  :      float
        @return                 :      returns avg. decoherence prob.
	"""
	param.Refresh()
	U = no.mixmatrix(param).U
	
	prob = 1.0/3.0
	#print (0.5*(U[a,0]**2-U[a,1]**2)*(U[b,0]**2-U[b,1]**2)+(1.0/6.0)*(U[a,0]**2+U[a,1]**2-2.0*U[a,2]**2)*(U[b,0]**2+U[b,1]**2-2.0*U[b,2]**2))
	prob = prob + np.exp(-2*alpha)*(0.5*(U[a,0]**2-U[a,1]**2)*(U[b,0]**2-U[b,1]**2)+(1.0/6.0)*(U[a,0]**2+U[a,1]**2-2.0*U[a,2]**2)*(U[b,0]**2+U[b,1]**2-2.0*U[b,2]**2))	
	return prob.real

#===============================================================================
# END DECOHERENCE FORMULAE
#===============================================================================

#===============================================================================
# BEGIN ADIABATIC FORMULAE
#===============================================================================

def Adiabaticity(x,E,dm2matter,sin2thmatter,Acc):
	# GIUNTI FUNDAMENTAL OF NEUTRINO PHYSICS p. 334
	# returns		: 	adiabaticity parameter
	# x				:	position
	# E				:	energy 
	# dm2matter		:	squared mass difference in matter
	# sin2thmatter	:	sin of 2(mixing angle) in matter
	# Acc			: 	function of the charge current potential
	dAcc = 0.0#sp.diff(a, n, axis)
	return (dm2matter**2)/(2*E*sin2thmatter*np.abs(dAcc(x)))

def SunAdiabaticity(dmsq,E):
	# Neutrino Physics Kai Suber, p. 279, eq. 10.89
	# E	: Neutrino Energy [eV]
	# dmsq : squared mass difference [eV^2]
	return 1.69e8*(dmsq/(E/pc.MeV))*pc.s12**2/pc.c12

def ParkeFormulaSurProb(E,track,body):
	# Neutrino Physics Kai Suber, p. 279, eq. 10.89
	# E	: Neutrino Energy [eV]
	# dmsq : squared mass difference [eV^2]
	Pc = np.exp(-pc.PIby2*SunAdiabaticity(pc.dm21sq, E))
	th = np.arcsin(pc.s12)/2.0
	s2m = sin2matter(pc.s12,pc.dm21sq,dm2matter(pc.s12,pc.dm21sq,Acc(E,track,body)))
	c2m = np.sqrt(1.0-s2m**2)
	#print Pc,np.sqrt(1.0-pc.s12**2)
	return np.sin(th)**2 + Pc*np.sqrt(1.0-pc.s12**2)

#===============================================================================
# END ADIABATIC FORMULAE
#===============================================================================


################################################################################################################
############################################ TWO NEUTRINO OSC PROB #############################################
################################################################################################################

def dm2matter(sin2thvacuum,dm2vacuum,Acc):
# GIUNTI FUNDAMENTAL OF NEUTRINO PHYSICS p. 333
# RETORNA DELTA DE MASA AL CUADRADO EN MATERIA 
# sin2thvacuum		: 	sin of 2(mixing angle) in vacuum
# dm2vacuum		:	squared mass diference in vacuum
# Acc			:	charge current matter potential
	cos2thvacuum = np.sqrt(1.0-sin2thvacuum**2)
	x = (dm2vacuum*cos2thvacuum-Acc)**2 + (dm2vacuum*sin2thvacuum)**2
	return np.sqrt(x)

def thmatter(thvacuum,dm2matter,Acc):
# GIUNTI FUNDAMENTAL OF NEUTRINO PHYSICS p. 333
# RETORNA ANGULO DE MEZCLA EN MATERIA
# Acc 		:	charge current matter potential  
# thvacuum	: 	mixing angle in vacuum
# dm2matter	: 	square matter difference in matter
	x = 1.0 - Acc/dm2matter/np.cos(2.0*thvacuum)
	tan2thm = np.tan(2.0*thvacuum)/x
	return np.arctan(tan2thm)/2.0

def sin2matter(sin2thvacuum,dm2vacuum,dm2matter):
# GIUNTI FUNDAMENTAL OF NEUTRINO PHYSICS p. 333
# returns 		:	sin(2theta23MATTER)
# sin2thvacuum	:	sin(2theta23vacuum)
# dm2vacuum 	:	squared mass diference in vacuum
# dm2matter 	: 	squared mass diference in matter
	return sin2thvacuum*(dm2vacuum/dm2matter)

def Acc_CC(param,E,body,track):
# charge current matter potential
# x		: position [km]
# E		: energy [eV]
# GF	: Fermi constant [cm^-3]
# assuming electron fraction as 0.5
	#print body.density(track)
	return 2.0*E*param.sqr2*param.GF*param.Na*param.cm**-3*body.density(track)*body.ye(track)	
	#return 2.0*np.sqrt(2.0)*E*pc.GF*pc.Na*pc.cm**-3*body.density(track)*0.7*0.59
	
def NeuOsc2g(a,b,E,track,body,param):
	""" Neutrino oscillation probability in 2-generation. Matter effects in constant density.
	
	Oscillation probability in two generation considering matter effects.
	#	E 	: neutrino energy 	[eV]
	#	L	: distance 			[eV^-1]
	#	body	: body 
	#	track	: trackbody
	#	a	: initial neutrino flavor (0 : e, 1 : mu, 2 : tau) # not really
	#	b	: final neutrino flavor (0 : e, 1 : mu, 2 : tau)   # not really
	"""
	param.Refresh()
	
	L = track.x
	track.x = L/2.0
	#print L
	#Acc = np.real(no.flavorAcc(param, E, body, track)[0,0])
	if param.neutype == "neutrino":
		Acc = Acc_CC(param,E,body,track)
	elif param.neutype == "antineutrino":
		Acc = -Acc_CC(param,E,body,track)
	else :
		print "Wrong neutype."
		quit()
	#print Acc
	#quit()
	
	dmM2 = dm2matter(np.sin(2.0*param.th[2,3]),param.dm2[2,3],Acc)
	#print dmM2, param.dm2[2,3]
	#thM  = thmatter(param.th[2,3],dmM2,Acc) # NOT WORKING MUST FIX!
	s2m = sin2matter(np.sin(2.0*param.th[2,3]),param.dm2[2,3],dmM2)
	s2m = np.sin(2.0*param.th[2,3])
	if a == b :
		return 1.0 - s2m**2*np.sin((dmM2*L)/(4*E))**2
		#return 1.0 - np.sin(2.0*thM)**2*np.sin((dmM2*L)/(4*E))**2
	else : 
		return 0.0 + s2m**2*np.sin((dmM2*L)/(4*E))**2
		#return 0.0 + np.sin(2.0*thM)**2*np.sin((dmM2*L)/(4*E))**2

################################################################################################################
########################################## THREE NEUTRINO OSC PROB #############################################
################################################################################################################

########################################### FORTRAN RK METHOD ##################################################

def neuosc3gfortran(enu,distance,sin2,dm,neu,matter):
	return float(subprocess.Popen('./probosc3gv3.exe '+ str(distance) + ' ' +  str(enu) + ' 0.0 0.0 2.9 ' + str(sin2)+ ' ' + str(dm) + ' ' + str(neu),shell=True,stdout=subprocess.PIPE).stdout.read())
	
def neuosc3gfortranNSI(enu,distance,sin2,dm,neu,matter,eps):
	return float(subprocess.Popen('./probnsi.exe '+ str(distance) + ' ' +  str(enu) + ' 0.0 0.0 2.9 ' + str(sin2)+ ' ' + str(dm) + ' ' + str(neu)+ ' ' + str(eps),shell=True,stdout=subprocess.PIPE).stdout.read())

def neuosc3(distance,enu,ow,ox,density,flavor):
	# flavor = 0 -> electron - electron
	# flavor = 1 -> electron - muon
	# flavor = 2 -> electron - tau
	try :
		return map(lambda x: float(x),subprocess.Popen('./probosc3g.exe '+ str(distance) + ' ' +  str(enu) + ' ' + str(ox) + ' ' + str(ow) + ' ' + str(density),shell=True,stdout=subprocess.PIPE).stdout.read().split())[flavor]
	except (ValueError,RuntimeError):
		# perturbative solution 
		eps = 0.001*np.random.uniform(0.0,1.0)
		#print eps
		return map(lambda x: float(x),subprocess.Popen('./probosc3g.exe '+ str(distance+eps) + ' ' +  str(enu-eps) + ' ' + str(ox+eps) + ' ' + str(ow-eps) + ' ' + str(density+eps),shell=True,stdout=subprocess.PIPE).stdout.read().split())[flavor]

########################################### ANALITIC CONSTANT DENSITY ##########################################

def NeuOsc3g(a,b,E,track,body,param):
	""" Formalae for standard 3-flavor neutrino oscillation for constant density case.
	# Using evolution operator formalism in constant density.
	#	E 	: neutrino energy 	[eV]
	#	track	: body track
	#	body	: body where neutrino propagates
	#	a	: initial neutrino flavor (0 : e, 1 : mu, 2 : tau)
	#	b	: final neutrino flavor (0 : e, 1 : mu, 2 : tau)
	"""
	fM2 = no.flavorM2(param)
	return no.OscProbConstantDensity(a,b,E,track,body,fM2,param)

def NeuOsc3gVacuum(a,b,E,L,param):
	""" Formalae for standard 3-flavor neutrino oscillation
	# FUNDAMENTALS OF NEUTRINO PHYSICS GIUNTI pag. 252
	#	E 	: neutrino energy 	[eV]
	#	L	: distance 		[eV^-1]
	#	a	: initial neutrino flavor (0 : e, 1 : mu, 2 : tau)
	#	b	: final neutrino flavor (0 : e, 1 : mu, 2 : tau)
	"""
	param.Refresh()
	PMNS = no.mixmatrix(param)
	U = PMNS.U
	UCT = PMNS.UCT
	rt,it =[],[]
	[[ rt.append(np.real(UCT[k,a]*U[b,k]*U[a,j]*UCT[j,b])*np.sin(param.dm2[j+1,k+1]*L/(4.0*E))**2)  for k in range(param.numneu) if k>j] for j in range(param.numneu)]
	[[ it.append(np.imag(UCT[k,a]*U[b,k]*U[a,j]*UCT[j,b])*np.sin(param.dm2[j+1,k+1]*L/(2.0*E)))     for k in range(param.numneu) if k>j] for j in range(param.numneu)]
	
	if a == b :
		return 1.0 - 4.0*np.sum(rt) + 2.0*np.sum(it)  
	else : 
		return 0.0 - 4.0*np.sum(rt) + 2.0*np.sum(it)
		
if __name__ == '__main__':
	pc = PC.PhysicsConstants()
	E = 115.0*pc.MeV
	L = 1500.0*pc.km
	import body as bd
	body = bd.Earth()
#	body = bd.Vacuum()
	track = body.track(L,L,L)
	pc.neutype = "antineutrino"
	print NeuOsc3g(0,0,E,track,body,pc)
#	print NeuOsc3gVacuum(0,0,E,L,pc)
