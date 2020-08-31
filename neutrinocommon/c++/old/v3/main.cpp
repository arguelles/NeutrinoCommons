#include <iostream>
#include "body.h"
#include "physconst.h"
#include "neuosc.h"

using namespace std;

int main()
{

PhysConst pc;

gsl_matrix_complex* fM2;
fM2 = flavorM2(pc);

Vacuum vacuum;
Vacuum::Track track(0.0,100.0*pc.km);

cout << "Track information" << endl;
cout << track.x << endl;
cout << track.xini << endl;
cout << track.xend << endl;

cout << "Body Information" << endl;
cout << vacuum.name << endl;
cout << vacuum.density(track) << endl;

double prob[pc.numneu];
double E = 10.0*pc.MeV;
CalNeuOscGSL(prob,0,E,track,vacuum,fM2,pc,1.0e-5,0.0,0,0);

cout << "Oscillation Probabilities" << endl;
for(int i =0;i<pc.numneu;i++){
    cout << prob[i] << endl;    
}
}
