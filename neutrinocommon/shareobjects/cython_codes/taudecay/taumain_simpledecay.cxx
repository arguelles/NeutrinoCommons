/**
 * Simulates tau decay and returns the energies of the
 * outgoing neutrinos.
 *
 * Modified from TAUOLA example. 
 *
 * @author C. Arguelles
 * @date 10 July 2011
 */

#include "taumain_simpledecay.h"

using namespace std;

HepMC::GenEvent * make_simple_tau_event(double tau_energy){
  /*
	IN  : 	tau_energy [GeV]
	OUT :	HepMCEvent 
	Creates an HepMC event that has a single tau of energy
	tau_energy.
  */
  //create HepMC event
  HepMC::GenEvent * event = new HepMC::GenEvent();
	
  //define tau mass
    double tau_mass = double(Tauolapp::parmas_.amtau);

  //make tau's four vector
  HepMC::FourVector momentum_tau1(0,0,0,0);

  momentum_tau1.setPz(sqrt(tau_energy*tau_energy-tau_mass*tau_mass));
  momentum_tau1.setE(tau_energy);

  //make particles
  HepMC::GenParticle * tau1 = new HepMC::GenParticle(momentum_tau1,15,1);
  
  //set the masses
  tau1->set_generated_mass(tau_mass);
 
  //make the vertex
  HepMC::GenVertex * vertex = new HepMC::GenVertex();
  vertex->add_particle_out(tau1);

  event->add_vertex(vertex);

  //calculate center of mass frame
  HepMC::FourVector cms(0,0,0,tau_mass);

  HepMC::GenParticle * cms_particle = new HepMC::GenParticle(cms,0,0);

  //Make TauolaParticles for doing boosting:
    Tauolapp::TauolaHepMCParticle * cms_boost = new Tauolapp::TauolaHepMCParticle(cms_particle);

    Tauolapp::TauolaHepMCParticle first_tau(tau1);

  first_tau.boostFromRestFrame(cms_boost);

  //clean up
  delete cms_boost;
  delete cms_particle;

  return event;
};

class IsNeutrino{
public:
	bool operator()(const HepMC::GenParticle* p){
		if (	(p->pdg_id() == 12 or p->pdg_id() == 14 or p->pdg_id() == 16 or p->pdg_id() == -12 or p->pdg_id() == -14 or p->pdg_id() == -16 ) 
			&& p->momentum().e() > 1.0		) return 1;
		return 0;
	}
};

class IsStateFinal{
public:
	bool operator()(const HepMC::GenParticle* p){
		if (!p->end_vertex() && p->status()==1) return 1;
		return 0;
	}
};

vector<double> NeutrinoEnergies(double E,int tau_init){

  int NumberOfEvents = 1;

  if (!tau_init){
	  //These three lines are not really necessary since they are the default
      Tauolapp::Tauola::setDecayingParticle(15);
//	  Tauolapp::Tauola::setEtaK0sPi(1,1,1);
      Tauolapp::Tauola::setSameParticleDecayMode(0);
      Tauolapp::Tauola::setOppositeParticleDecayMode(0);
      Tauolapp::Tauola::initialize();
  }

  // Begin event loop. Generate event.
  for (int iEvent = 0; iEvent < NumberOfEvents; ++iEvent) {

	// Convert event record to TAUOLA format
//        double energy = 0.0;
	HepMC::GenEvent * event = make_simple_tau_event(E);
      Tauolapp::TauolaHepMCEvent * t_event = new Tauolapp::TauolaHepMCEvent(event);

	// Decay taus
	t_event->decayTaus();
	// print event
	event->print();
	
	// get final state particles that are neutrinos
	HepMC::GenEvent* evt = t_event->getEvent();

//	IsStateFinal isfinal;
//	IsNeutrino isneutrino;

	vector<double> neu_type_energy;

	std::list<HepMC::GenParticle*> end_neutrinos;
	for ( HepMC::GenEvent::particle_iterator p = evt->particles_begin();
            p != evt->particles_end(); ++p ) {
//		cout << (*p)->pdg_id() << " " <<(*p)->momentum().e()<<endl;
		if( (*p)->pdg_id() == 16 or (*p)->pdg_id() == -16 or(*p)->pdg_id() == 14 or(*p)->pdg_id() == -14 or(*p)->pdg_id() == 12 or(*p)->pdg_id() == -12){
//		if( isneutrino(*p) ){
//		energy = (*p)->momentum().e();
		neu_type_energy.push_back(double((*p)->pdg_id()));
		neu_type_energy.push_back(double((*p)->momentum().e()));
//		cout << (*p)->pdg_id() << " " <<(*p)->momentum().e()<<endl;
		};
	}

    //clean up
    delete event;
    delete t_event;

//    return energy;
    return neu_type_energy;
  }
};

