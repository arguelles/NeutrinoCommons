#ifndef __MLINTERFACE_H
#define __MLINTERFACE_H

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