#include <iostream>
#include "body.h"
#include "physconst.h"
#include "neuosc.h"

using namespace std;

int main()
{

Vacuum vacuum("hola",5.0);
Vacuum::Track track(0.0,1.0);

cout << track.x << endl;
cout << track.xini << endl;
cout << track.xend << endl;

cout << vacuum.name << endl;
cout << vacuum.density(track) << endl;

PhysConst pc;
cout << pc.keV << endl;

cout << gsl_matrix_get(pc.th,1,2) << endl;
cout << gsl_matrix_get(pc.c,1,2) << endl;
pc.th12 = 0.36;
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
}
