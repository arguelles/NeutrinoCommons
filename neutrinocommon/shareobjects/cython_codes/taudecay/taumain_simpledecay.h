#ifndef _tausimpledecay_h_included_
#define _tausimpledecay_h_included_

//C++ std headers
#include <list>
#include <vector>

//TAUOLA headers
#include "Tauola.h"
#include "TauolaHepMCEvent.h"

//HepMC headers
#include "HepMC/GenEvent.h"
#include "HepMC/GenParticle.h"


HepMC::GenEvent * make_simple_tau_event(double);

vector<double> NeutrinoEnergies(double,int);

#endif  
