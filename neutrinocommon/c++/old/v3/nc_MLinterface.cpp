#include "nc_MLinterface.h"

//----------------------------------------------------------------------
//------------------------ MATHEMATICA INTERFACE -----------------------
//----------------------------------------------------------------------

//Definitions for argp
//const char *argp_program_version = "NeutrinoCommon-C++ MLInterface 0.1";
//const char *argp_program_bug_address = "<carguellesdel@gmail.com>";
//static char argp_doc[] = "Neutrino oscillation probability code for Mathematica.";
//static char argp_option_doc[] = "[options]";

int main(int argc, char **argv){
    int status = -1;
    // setting defaults parameters
    SetDefaultParameters();
    
    // command line arguments
    // Parse command line arguments
    //struct argp argp = { cmdline_options, parse_opt, argp_option_doc, argp_doc };
    //argp_parse(&argp, argc, argv, ARGP_SILENT, NULL, NULL);
    
    // Run in Mathlink mode
    status = MLMain(argc,argv);
    
    return status;
}

//-----------------------------------------------------------------------
//FUNCTIONS
//-----------------------------------------------------------------------

int SetDefaultParameters(void){
    return 1;
}

int SetParams(string param_name,double param_value){
    return 1;
}

void GetParams(void){
    
}

void CalNeuOscGSLML(int ineu,double E,string body_name,double body_params[],long body_params_len,double track_params[],long track_params_len,string neutype,double abs_error,double rel_error,int return_state,int optimization){
    int numneu = GLB_param.numneu;
    double osc_prob[numneu];
    Body body;
    Track track;
    
    if (body_name == "Vacuum"){
        // checking parameters lengths
        if (body_params_len == 0){
            Body body;    
        } else {
            MLEvaluate(stdlink, "Message[CalNeuOscGSLML::Invalid body arguments.]");
        }
        
        if (track_params_len == 2){
            Track trackInit(track_params[0],track_params[1]);
            track = trackInit;
        } else {
            MLEvaluate(stdlink, "Message[CalNeuOscGSLML::Invalid track arguments.]");
        }
    } else {
        MLEvaluate(stdlink, "Message[CalNeuOscGSLML::Invalid body name.]");
        exit(0);
    }
    
    PhysConst param = GLB_param;
    gsl_matrix_complex* fM2 = flavorM2(GLB_param);

    CalNeuOscGSL(osc_prob,ineu,E,track,body,fM2,param,abs_error,rel_error,return_state,optimization);
    
    // passing the probability back to Mathematica
    if (return_state){
        MLPutRealList(stdlink, osc_prob, 2*numneu);
    } else {
        MLPutRealList(stdlink, osc_prob, numneu);
    }
}