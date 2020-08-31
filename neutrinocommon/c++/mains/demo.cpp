/*
Author  : C.A.Arguelles
Date    : 3/JUL/2012

This files explains the basic usage of the neutrinocommon C++
Runge Kutta.

*/

#include <iostream>
#include "body.h"
#include "physconst.h"
#include "neuosc.h"

using namespace std;

int main()
{

/*
  This creates the analogus physicsconstants
  object.
*/

PhysConst pc;

/*
  We can calculate the square mass difference
  in the flavor base corresponding to this
  parameter set.
*/

gsl_matrix_complex* fM2;
fM2 = flavorM2(pc);

/*
  We must define the body on which the neutrino
  propagates and its corresponding trayectory.
*/

/*
 For vaccum
*/

Vacuum vacuum;
Vacuum::Track vacuum_track(0.0,100.0*pc.km);

cout << "Track information" << endl;
cout << vacuum_track.x << endl;
cout << vacuum_track.xini << endl;
cout << vacuum_track.xend << endl;

cout << "Body Information" << endl;
cout << vacuum.name << endl;
cout << vacuum.density(vacuum_track) << endl;

/*
  We have to create an array where the
  calculation will be stored.
*/

double prob[pc.numneu];

/*
  Setting the neutrino energy in eV.
*/

double E = 100.0*pc.MeV;

/*
  We calculate the oscillarion probability.
    CalNeuOscGSL(probability,
                initial_neutrino_flavor,
                neutrino_energy,
                neutrino_trayectory,
                body,
                square_mass_matrix in flavor basis,
                physics parameters,
                absolute error,
                relative error,
                return_neutrino_state,
                use optimizations)
*/
CalNeuOscGSL(prob,0,E,vacuum_track,vacuum,fM2,pc,1.0e-5,1.0e-5,0,0);

/*
  Printing the oscillation probabilities.
*/

cout << "Vacuum Oscillation Probabilities" << endl;
for(int i =0;i<pc.numneu;i++){
    cout << prob[i] << endl;    
}

/*
  We can make similar calculations for other
  bodies.
*/

/* Earth */
double baseline = 500.0*pc.km;

Earth earth;
Earth::Track earth_track(0.0,baseline/2.0,baseline);

CalNeuOscGSL(prob,0,E,earth_track,earth,fM2,pc,1.0e-5,0.0,0,0);

cout << "Earth Oscillation Probabilities" << endl;
for(int i =0;i<pc.numneu;i++){
    cout << prob[i] << endl;    
}

/* Constant Density */

ConstantDensity constdens(10.0,0.5);
ConstantDensity::Track constdens_track(0.0,baseline);

CalNeuOscGSL(prob,0,E,constdens_track,constdens,fM2,pc,1.0e-5,0.0,0,0);

cout << "Constant Density Oscillation Probabilities" << endl;
for(int i =0;i<pc.numneu;i++){
    cout << prob[i] << endl;    
}

/* Variable Density */

VariableDensity vardens(10.0,2.0,0.5);
VariableDensity::Track vardens_track(0.0,baseline);

CalNeuOscGSL(prob,0,E,vardens_track,vardens,fM2,pc,1.0e-5,0.0,0,0);

cout << "Variable Density Oscillation Probabilities" << endl;
for(int i =0;i<pc.numneu;i++){
    cout << prob[i] << endl;    
}

/* Star Density */
double star_radius = 1.0e5*pc.km;
double ini_radius = 1000.0*pc.km;
double ini_density = 100.0*pc.km;

Star star(star_radius,ini_radius,ini_density);
Star::Track star_track(0.0,star_radius);

CalNeuOscGSL(prob,0,E,star_track,star,fM2,pc,1.0e-5,0.0,0,0);

cout << "Star Oscillation Probabilities" << endl;
for(int i =0;i<pc.numneu;i++){
    cout << prob[i] << endl;    
}
}
