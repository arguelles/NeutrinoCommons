#ifndef __PHYSCONST_H
#define __PHYSCONST_H

#include <string>
#include <gsl/gsl_matrix.h>
#include <math.h>
using namespace std;

class PhysConst {
    public : 
        // class identifiers
        string name; 
        string linestyle;
        string markerstyle;
        string colorstyle;
        string savefilename;
        // mathematical constants //
        double pi;
        double piby2; 
        double sqrt2;
        double ln2;
        // astronomical constants //
        double earthradius;
        double sunradius;
        ///// physics constants/////
        double GF;
        double Na;
        double sw_sq;
        double G;
        double alpha;
        /////////// units //////////
        // energy
        double TeV;
        double GeV;
        double MeV;
        double keV;
        double Joule;
        // mass
        double kg;
        double gr;
        // time 
        double sec;
        double hour;
        double day;
        double year;
        // distance
        double meter;
        double cm;
        double km;
        double fermi;
        double angstrom;
        double AU;
        double parsec;
        // luminocity
        double picobarn;
        double femtobarn;
        // presure
        double Pascal;
        double hPascal;
        double atm;
        double psi;
        // temperature
        double Kelvin;
        // angle
        double degree;
        ////// neutrino osc. param. //////////
        // basic
        int numneu; 
        int numneumax;
        string neutype;
        // mixing angles
        double th12;
        double th13;
        double th23;
        double th14;
        double th24;
        double th34;
        double th15;
        double th25;
        double th35;
        double th45;
        double th16;
        double th26;
        double th36;
        double th46;
        double th56;
        // square mass differences
        double dm21sq;
        double dm31sq;
        double dm41sq;
        double dm51sq;
        double dm61sq;
        // cp-phases
        double delta1;
        double delta2;
        double delta3;
        // matrices
        // angles
        gsl_matrix *th;
        // cosine mixing matrix
        gsl_matrix *c;
        // sine mixing matrix
        gsl_matrix *s;
        // cp-phases
        gsl_matrix *dcp;
        // square mass differences
        gsl_matrix *dmsq;
        // mixing matrix
        gsl_matrix_complex *U;
        
        int Refresh(void);
        PhysConst(void);
};

#endif
