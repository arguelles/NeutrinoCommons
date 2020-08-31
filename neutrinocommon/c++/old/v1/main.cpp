#include <iostream>
#include "body.h"
#include "physconst.h"
#include "neuosc.h"

using namespace std;

int main()
{

PhysConst pc;
cout << pc.keV << endl;

cout << gsl_matrix_get(pc.th,1,2) << endl;
cout << gsl_matrix_get(pc.c,1,2) << endl;
//pc.th12 = 0.36;
pc.Refresh();
cout << gsl_matrix_get(pc.c,1,2) << endl;

gsl_matrix_complex* rotation;
rotation = R(1,2,0,pc);

cout << "R_12 rotation matrix" << endl;
gsl_matrix_complex_fprintf(stdout,rotation,"%f");

gsl_matrix_complex* mixmatrix;
mixmatrix = MixMatrix(pc);

cout << "Mixing matrix" << endl;
gsl_matrix_complex_fprintf(stdout,mixmatrix,"%f");

gsl_matrix_complex* fM2;
fM2 = flavorM2(pc);
cout << "Square mass diff. in flavor basis" << endl;
gsl_matrix_complex_fprintf(stdout,fM2,"%f");

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
double E = 10.0*pc.GeV;
CalNeuOscGSL(prob,0,E,track,vacuum,fM2,pc,1.0e-5,0.0,0,0);

cout << "Oscillation Probabilities" << endl;
for(int i =0;i<pc.numneu;i++){
    cout << prob[i] << endl;    
}
}
