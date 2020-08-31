#include "physconst.h"

PhysConst::PhysConst(void){
    /* PHYSICS CONSTANTS
    #===========================================================================
    # NAME
    #===========================================================================
    */
    
    name = "STD";                    // Default values
    linestyle = "solid";             // Default linestyle in plots
    markerstyle = "*";               // Default marker style
    colorstyle = "red";              // Default color style
    savefilename = "output.dat";     // Default color style
    
    /*
    #===============================================================================
    # MATH
    #===============================================================================
    */
    pi=3.14159265;		    // Pi
    piby2=1.5707963268;	            // Pi/2
    sqrt2=1.4142135624;	            // Sqrt[2]
    ln2 = log(2.0);                 // log[2]
    
    /*
    #===============================================================================
    # EARTH 
    #===============================================================================
    */
    earthradius = 6371.0;	    // [km] Earth radius
    /*
    #===============================================================================
    # SUN 
    #===============================================================================
    */
    sunradius = 109.0*earthradius;  // [km] Sun radius 
    
    /*
    #===============================================================================
    # # PHYSICAL CONSTANTS
    #===============================================================================
    */
    GF = 1.16639e-23;	            // [eV^-2] Fermi Constant 
    Na = 6.0221415e+23;		    // [mol cm^-3] Avogadro Number
    sw_sq = 0.2312;                 // [dimensionless] sin(th_weinberg) ^2
    G  = 6.67300e-11;               // [m^3 kg^-1 s^-2]
    alpha = 1.0/137.0;              // [dimensionless] fine-structure constant 
    
    /*
    #===============================================================================
    # UNIT CONVERSION FACTORS
    #===============================================================================
    */
    // Energy
    TeV = 1.0e12;                   // [eV/TeV]
    GeV = 1.0e9;                    // [eV/GeV]
    MeV = 1.0e6;                    // [eV/MeV]
    keV = 1.0e3;                    // [eV/keV]
    Joule = 1/1.60225e-19;          // [eV/J]
    // Mass
    kg = 5.62e35;                   // [eV/kg]
    gr = 1e-3*kg;                   // [eV/g] 
    // Time
    sec = 1.523e15;                 // [eV^-1/s]
    hour = 3600.0*sec;              // [eV^-1/h]
    day = 24.0*hour;                // [eV^-1/d]
    year = 365.0*day;               // [eV^-1/yr]
    // Distance
    meter = 5.076e6;                // [eV^-1/m]
    cm = 1.0e-2*meter;              // [eV^-1/cm]
    km = 1.0e3*meter;               // [eV^-1/km]
    fermi = 1.0e-15*meter;          // [eV^-1/fm]
    angstrom = 1.0e-10*meter;       // [eV^-1/A]
    AU = 149.60e9*meter;            // [eV^-1/AU]
    parsec = 3.08568025e16*meter;   // [eV^-1/parsec]
    // luminocity
    picobarn = 1.0e-36*pow(cm,2);       // [eV^-2/pb]
    femtobarn = 1.0e-39*pow(cm,2);      // [eV^-2/fb]
    // Presure
    Pascal = Joule/pow(meter,3);        // [eV^4/Pa]
    hPascal = 100.0*Pascal;         // [eV^4/hPa]
    atm = 101325.0*Pascal;          // [eV^4/atm]
    psi = 6893.0*Pascal;            // [eV^4/psi]
    // Temperature
    Kelvin = 1/1.1604505e4;         // [eV/K]
    // Angle
    degree = pi/180.0;              // [rad/degree]
    
    /*
    #===============================================================================
    # NEUTRINO OSCILLATION PARAMETERS 
    #===============================================================================
    */
    
    numneu = 3;                     // number of neutrinos
    numneumax = 6;                  // maximum neutrino number
    neutype = "neutrino";           // neutrino or antineutrino
    
    // angles
    th12 = 0.563942;
    th13 = 0.154085;
    th23 = piby2/2.0;
    th14 = 0.0;
    th24 = 0.0;
    th34 = 0.0;
    th15 = 0.0;
    th25 = 0.0;
    th35 = 0.0;
    th45 = 0.0;
    th16 = 0.0;
    th26 = 0.0;
    th36 = 0.0;
    th46 = 0.0;
    th56 = 0.0;
    
    // square mass differences
    dm21sq = 7.65e-5;
    dm31sq = 2.47e-3;
    dm41sq = 0.0;
    dm51sq = 0.0;
    dm61sq = 0.0;
    
    // cp-phases
    delta1 = 0.0;
    delta2 = 0.0;
    delta3 = 0.0;
    
    // initializing matrices
    
    dmsq = gsl_matrix_alloc(numneumax,1);
    gsl_matrix_set(dmsq,0,0,0.0);
    gsl_matrix_set(dmsq,1,0,dm21sq);
    gsl_matrix_set(dmsq,2,0,dm31sq);
    gsl_matrix_set(dmsq,3,0,dm41sq);
    gsl_matrix_set(dmsq,4,0,dm51sq);
    gsl_matrix_set(dmsq,5,0,dm61sq);
    
    th = gsl_matrix_alloc(numneumax+1,numneumax+1);
    gsl_matrix_set(th,1,2,th12);
    gsl_matrix_set(th,1,2,th12);
    gsl_matrix_set(th,1,3,th13);
    gsl_matrix_set(th,2,3,th23);
    gsl_matrix_set(th,1,4,th14);
    gsl_matrix_set(th,2,4,th24);
    gsl_matrix_set(th,3,4,th34);
    gsl_matrix_set(th,1,5,th15);
    gsl_matrix_set(th,2,5,th25);
    gsl_matrix_set(th,3,5,th35);
    gsl_matrix_set(th,4,5,th45);
    gsl_matrix_set(th,1,6,th16);
    gsl_matrix_set(th,2,6,th26);
    gsl_matrix_set(th,3,6,th36);
    gsl_matrix_set(th,4,6,th46);
    gsl_matrix_set(th,5,6,th56);
    
    c = gsl_matrix_alloc(numneumax+1,numneumax+1);
    gsl_matrix_set(c,1,2,cos(th12));
    gsl_matrix_set(c,1,3,cos(th13));
    gsl_matrix_set(c,1,4,cos(th14));
    gsl_matrix_set(c,2,3,cos(th23));
    gsl_matrix_set(c,2,4,cos(th24));
    gsl_matrix_set(c,3,4,cos(th34));
    gsl_matrix_set(c,1,5,cos(th15));
    gsl_matrix_set(c,2,5,cos(th25));
    gsl_matrix_set(c,3,5,cos(th35));
    gsl_matrix_set(c,4,5,cos(th45));
    gsl_matrix_set(c,1,6,cos(th16));
    gsl_matrix_set(c,2,6,cos(th26));
    gsl_matrix_set(c,3,6,cos(th36));
    gsl_matrix_set(c,4,6,cos(th46));
    gsl_matrix_set(c,5,6,cos(th56));
    
    s = gsl_matrix_alloc(numneumax+1,numneumax+1);
    gsl_matrix_set(s,1,2,sin(th12));
    gsl_matrix_set(s,1,3,sin(th13));
    gsl_matrix_set(s,1,4,sin(th14));
    gsl_matrix_set(s,2,3,sin(th23));
    gsl_matrix_set(s,2,4,sin(th24));
    gsl_matrix_set(s,3,4,sin(th34));
    gsl_matrix_set(s,1,5,sin(th15));
    gsl_matrix_set(s,2,5,sin(th25));
    gsl_matrix_set(s,3,5,sin(th35));
    gsl_matrix_set(s,4,5,sin(th45));
    gsl_matrix_set(s,1,6,sin(th16));
    gsl_matrix_set(s,2,6,sin(th26));
    gsl_matrix_set(s,3,6,sin(th36));
    gsl_matrix_set(s,4,6,sin(th46));
    gsl_matrix_set(s,5,6,sin(th56));      
    
    dcp = gsl_matrix_alloc(numneumax-2+1,1);
    gsl_matrix_set(dcp,0,0,1.0);
    gsl_matrix_set(dcp,1,0,delta1);
    gsl_matrix_set(dcp,2,0,delta2);
    gsl_matrix_set(dcp,3,0,delta3);        
};

int PhysConst::Refresh(void){
    // reinitializing matrices
    
    dmsq = gsl_matrix_alloc(numneumax,1);
    gsl_matrix_set(dmsq,0,0,0.0);
    gsl_matrix_set(dmsq,1,0,dm21sq);
    gsl_matrix_set(dmsq,2,0,dm31sq);
    gsl_matrix_set(dmsq,3,0,dm41sq);
    gsl_matrix_set(dmsq,4,0,dm51sq);
    gsl_matrix_set(dmsq,5,0,dm61sq);
    
    th = gsl_matrix_alloc(numneumax+1,numneumax+1);
    gsl_matrix_set(th,1,2,th12);
    gsl_matrix_set(th,1,2,th12);
    gsl_matrix_set(th,1,3,th13);
    gsl_matrix_set(th,2,3,th23);
    gsl_matrix_set(th,1,4,th14);
    gsl_matrix_set(th,2,4,th24);
    gsl_matrix_set(th,3,4,th34);
    gsl_matrix_set(th,1,5,th15);
    gsl_matrix_set(th,2,5,th25);
    gsl_matrix_set(th,3,5,th35);
    gsl_matrix_set(th,4,5,th45);
    gsl_matrix_set(th,1,6,th16);
    gsl_matrix_set(th,2,6,th26);
    gsl_matrix_set(th,3,6,th36);
    gsl_matrix_set(th,4,6,th46);
    gsl_matrix_set(th,5,6,th56);
    
    c = gsl_matrix_alloc(numneumax+1,numneumax+1);
    gsl_matrix_set(c,1,2,cos(th12));
    gsl_matrix_set(c,1,3,cos(th13));
    gsl_matrix_set(c,1,4,cos(th14));
    gsl_matrix_set(c,2,3,cos(th23));
    gsl_matrix_set(c,2,4,cos(th24));
    gsl_matrix_set(c,3,4,cos(th34));
    gsl_matrix_set(c,1,5,cos(th15));
    gsl_matrix_set(c,2,5,cos(th25));
    gsl_matrix_set(c,3,5,cos(th35));
    gsl_matrix_set(c,4,5,cos(th45));
    gsl_matrix_set(c,1,6,cos(th16));
    gsl_matrix_set(c,2,6,cos(th26));
    gsl_matrix_set(c,3,6,cos(th36));
    gsl_matrix_set(c,4,6,cos(th46));
    gsl_matrix_set(c,5,6,cos(th56));
    
    s = gsl_matrix_alloc(numneumax+1,numneumax+1);
    gsl_matrix_set(s,1,2,sin(th12));
    gsl_matrix_set(s,1,3,sin(th13));
    gsl_matrix_set(s,1,4,sin(th14));
    gsl_matrix_set(s,2,3,sin(th23));
    gsl_matrix_set(s,2,4,sin(th24));
    gsl_matrix_set(s,3,4,sin(th34));
    gsl_matrix_set(s,1,5,sin(th15));
    gsl_matrix_set(s,2,5,sin(th25));
    gsl_matrix_set(s,3,5,sin(th35));
    gsl_matrix_set(s,4,5,sin(th45));
    gsl_matrix_set(s,1,6,sin(th16));
    gsl_matrix_set(s,2,6,sin(th26));
    gsl_matrix_set(s,3,6,sin(th36));
    gsl_matrix_set(s,4,6,sin(th46));
    gsl_matrix_set(s,5,6,sin(th56));      
    
    dcp = gsl_matrix_alloc(numneumax-2+1,1);
    gsl_matrix_set(dcp,0,0,1.0);
    gsl_matrix_set(dcp,1,0,delta1);
    gsl_matrix_set(dcp,2,0,delta2);
    gsl_matrix_set(dcp,3,0,delta3);    
    
    return 0;
}