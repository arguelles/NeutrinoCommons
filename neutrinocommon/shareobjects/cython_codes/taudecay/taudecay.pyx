from libcpp.vector cimport vector

cdef extern from "taumain_simpledecay.cxx":
#cdef extern from "taumain_simpledecay.h":
	vector[double] NeutrinoEnergies(double,int)

cpdef TaudecayNeutrinos(E,tauolainitialize = 0):
	cdef unsigned int i
	neu_type_energy = []
	neu_vec = NeutrinoEnergies(E,tauolainitialize)
	for i in range(neu_vec.size()):
		neu_type_energy.append(neu_vec[i])
	return neu_type_energy
