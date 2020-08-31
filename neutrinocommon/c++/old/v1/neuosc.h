#ifndef __NEUOSC_H
#define __NEUOSC_H

#include <iostream>
#include "body.h"
#include "physconst.h"
#include <float.h>
#include <math.h>
#include <complex.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_odeiv.h>

//#define CalNeuOscGSL_DEBUG
//#define RHS_GSL_DEBUG
//#define flavorH_DEBUG
//#define flavorM2_DEBUG
//#define flavorAcc_DEBUG
//#define MixMatrix_DEBUG

// gsl-rk auxiliary class
typedef struct
{
    public :
    // basic storage
    Track track;
    PhysConst param;
    gsl_matrix_complex *flavorM2;
    Body body;
    double E;
    
    // used for optimization
    gsl_matrix_complex H0;
    gsl_matrix_complex S;
    double rpos;
    gsl_matrix_complex cH;
    gsl_matrix_complex cS;
} Container;

// mixing and rotation matrices
gsl_matrix_complex* R(int,int,int,PhysConst);
gsl_matrix_complex* MixMatrix(PhysConst);

// hamiltonian definition
gsl_matrix_complex* massM2(PhysConst);
gsl_matrix_complex* flavorM2(PhysConst);
gsl_matrix_complex* flavorAcc(PhysConst,double,Body,Track);
gsl_matrix_complex* flavorH(PhysConst,double,Body,Track,gsl_matrix_complex*);

// rhs
int RHS_INT_GSL(double,gsl_vector,void *);
static int RHS_GSL(double,const double *, double * ,void *);

// runge kutta implementation
int CalNeuOscGSL(double [],int,double,Track,Body,gsl_matrix_complex*,PhysConst,double,double,bool,bool);
#endif
