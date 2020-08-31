#ifndef __MLINTERFACE_H
#define __MLINTERFACE_H

#include "physconst.h"
#include "body.h"
#include "global.h"
#include "neuosc.h"
#include <mathlink.h>
#include <argp.h>

//----------------------------------------------------------------------
//------------------------ MATHEMATICA INTERFACE -----------------------
//----------------------------------------------------------------------

// parameters
int SetDefaultParameters(void);
int SetParams(string,double);
void GetParams(void);
// oscillation Calculation
void CalNeuOscGSLML(int,double,string,double [],long,double [],long,string,double,double,int,int)

#endif