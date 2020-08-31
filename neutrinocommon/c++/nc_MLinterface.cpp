//----------------------------------------------------------------------
//------------------------ MATHEMATICA INTERFACE -----------------------
//----------------------------------------------------------------------

//Definitions for argp
//const char *argp_program_version = "NeutrinoCommon-C++ MLInterface 0.1";
//const char *argp_program_bug_address = "<carguellesdel@gmail.com>";
//static char argp_doc[] = "Neutrino oscillation probability code for Mathematica.";
//static char argp_option_doc[] = "[options]";

#define __MLDEBUG

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <argp.h>
#ifndef __MLDEBUG
#include <mathlink.h>
#endif
#include "physconst.h"
#include "body.h"
#include "global.h"
#include "neuosc.h"

//extern PhysConst GLB_param;

//-----------------------------------------------------------------------
//FUNCTIONS
//-----------------------------------------------------------------------

int SetDefaultParameters(){
    return 1;
}

int SetParams(int param_id,double param_value){
    return 1;
}

void GetParams(){
    
}

void CalNeuOscGSLML(int ineu,double E,int body_id,double body_params[],long int body_params_len,double track_params[],long int track_params_len,int neutype,double abs_error,double rel_error,int return_state,int optimization){
    int numneu = GLB_param.numneu;
    double osc_prob[numneu];
    Body body;
    Track track;
    
    if (body_id == 0){
    //body_id = 0 : Vacuum
        // checking parameters lengths
        if (body_params_len == 0){
            Body body;    
        } else {
            #ifndef __MLDEBUG
            MLEvaluate(stdlink, "Message[CalNeuOscGSLML::Invalid body arguments.]");
            #endif
        }
        
        if (track_params_len == 2){
            Track trackInit(track_params[0],track_params[1]);
            track = trackInit;
        } else {
            #ifndef __MLDEBUG
            MLEvaluate(stdlink, "Message[CalNeuOscGSLML::Invalid track arguments.]");
            #endif
        }
    } else {
        #ifndef __MLDEBUG
        MLEvaluate(stdlink, "Message[CalNeuOscGSLML::Invalid body name.]");
        #endif
        exit(0);
    }
    
    PhysConst param = GLB_param;
    gsl_matrix_complex* fM2 = flavorM2(GLB_param);
    
    if (neutype == 0){
        GLB_param.neutype == "neutrino";
    } else {
        GLB_param.neutype == "antineutrino";
    }

    CalNeuOscGSL(osc_prob,ineu,E,track,body,fM2,param,abs_error,rel_error,return_state,optimization);
    
    // passing the probability back to Mathematica
    if (return_state){
        #ifndef __MLDEBUG
        MLPutRealList(stdlink, osc_prob, 2*numneu);
        #endif
    } else {
        #ifndef __MLDEBUG
        MLPutRealList(stdlink, osc_prob, numneu);
        #endif
    }
    
    #ifdef __MLDEBUG
    printf("%g %g %g \n",osc_prob[0],osc_prob[1],osc_prob[2]);
    #endif
}

int main(int argc, char **argv){
    int status = -1;
    // setting defaults parameters
    SetDefaultParameters();
    
    // command line arguments
    // Parse command line arguments
    //struct argp argp = { cmdline_options, parse_opt, argp_option_doc, argp_doc };
    //argp_parse(&argp, argc, argv, ARGP_SILENT, NULL, NULL);
    
    // Run in Mathlink mode
    #ifndef __MLDEBUG
    status = MLMain(argc,argv);
    #endif
    
    #ifdef __MLDEBUG
    int ineu = 0;
    double E = 1000.0;
    int body_id = 0;
    double body_params[0];
    body_params[0] = 0.0;
    int body_params_len = 1;
    double track_params[2];
    track_params[0] = 0.0;
    track_params[1] = 10.0;
    int track_params_len = 2;
    int neutype = 0;
    double abs_error = 1.0e-5;
    double rel_error = 0.0;
    int return_state = 0;
    int optimization = 0;
    CalNeuOscGSLML(ineu,E,body_id,body_params,body_params_len,track_params,track_params_len,neutype,abs_error,rel_error,return_state,optimization);
    #endif   
    return status;
}

