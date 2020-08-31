// -------------------------------------------
// neutrinocommon.tm : MathLink template file
// -------------------------------------------

// Function declarations

// Physics Parameters Management

:Begin:
:Function:      SetDefaultParameters
:Pattern:       SetDefaultParameters[]
:Arguments:     {}
:ArgumentTypes: {}
:ReturnType:    Integer
:End:

:Begin:
:Function:      SetParamsML
:Pattern:       SetParamsML[name_Integer,value_?NumericQ]
:Arguments:     {name,value}
:ArgumentTypes: {Integer,Real}
:ReturnType:    Integer
:End:

:Begin:
:Function:      GetParamsML
:Pattern:       GetParamsML[]
:Arguments:     {}
:ArgumentTypes: {}
:ReturnType:    Manual
:End:

// Oscillation Calculation

:Begin:
:Function:      CalNeuOscGSLML
:Pattern:       CalNeuOscGSLML[ineu_Integer,energy_?NumericQ,bodyname_Integer,bodyparams_List,trackparams_List,neutype_Integer,abserror_?NumericQ,relerror_?NumericQ,returnstate_Integer,optimization_Integer]
:Arguments:     {ineu,energy,bodyname,bodyparams,trackparams,neutype,abserror,relerror,returnstate,optimization}
:ArgumentTypes: {Integer,Real,Integer,RealList,RealList,Integer,Real,Real,Integer,Integer}
:ReturnType:    Manual
:End:


#include <stdio.h>
#include <stdlib.h>
#include <mathlink.h>
#include <math.h>
#include <argp.h>
#include "physconst.h"
#include "body.h"
#include "neuosc.h"

PhysConst GLB_param;

//-----------------------------------------------------------------------
//FUNCTIONS
//-----------------------------------------------------------------------

int SetDefaultParameters(){

    GLB_param.numneu = 3;
    GLB_param.th12 = 0.563943;
    GLB_param.th13 = 0.154085;
    GLB_param.th23 = 1.570796;
    GLB_param.dm21sq = 7.65e-5;
    GLB_param.dm31sq = 2.47e-3;
    GLB_param.delta1 = 1.570796;
    
    GLB_param.Refresh();
    return 1;
}

int SetParamsML(int param_id,double param_value){
    // defining parameter id order
    switch (param_id)
    {
    case 0 :
        GLB_param.numneu = (int) param_value;
        break;
    case 1 :
        GLB_param.th12 = param_value;
        break;
    case 2 :
        GLB_param.th13 = param_value;
        break;
    case 3 :
        GLB_param.th23 = param_value;
        break;
    case 4 :
        GLB_param.th14 = param_value;
        break;
    case 5 :
        GLB_param.th24 = param_value;
        break;
    case 6 : 
        GLB_param.th34 = param_value;
        break;
    case 7 :
        GLB_param.th15 = param_value;
        break;
    case 8 :
        GLB_param.th25 = param_value;
        break;
    case 9 :
        GLB_param.th35 = param_value;
        break;
    case 10 :
        GLB_param.th45 = param_value;
        break;
    case 11 :
        GLB_param.th16 = param_value;
        break;
    case 12 :
        GLB_param.th26 = param_value;
        break;
    case 13 :
        GLB_param.th36 = param_value;
        break;
    case 14 :
        GLB_param.th46 = param_value;
        break;
    case 15 :
        GLB_param.th56 = param_value;
        break;
    case 16 :
        GLB_param.dm21sq = param_value;
        break;
    case 17 :
        GLB_param.dm31sq = param_value;
        break;
    case 18 :
        GLB_param.dm41sq = param_value;
        break;
    case 19 :
        GLB_param.dm51sq = param_value;
        break;
    case 20 :
        GLB_param.dm61sq = param_value;
        break;
    case 21 :
        GLB_param.delta1 = param_value;
        break;
    case 22 :
        GLB_param.delta2 = param_value;
        break;
    case 23 :
        GLB_param.delta3 = param_value;
        break;
        
    default:
            MLEvaluate(stdlink, "Message[SetParams::Invalid parameter id.]");
            MLClearError(stdlink);
            MLNextPacket(stdlink);
            MLNewPacket(stdlink);
            MLPutSymbol(stdlink, "$Failed");
    };
    
    GLB_param.Refresh();
    return 1;
}

void GetParamsML(){
    int paramnum = 24;
    double parameters[paramnum];

    parameters[0] = (double) GLB_param.numneu;
    parameters[1] = GLB_param.th12;
    parameters[2] = GLB_param.th13;
    parameters[3] = GLB_param.th23;
    parameters[4] = GLB_param.th14;
    parameters[5] = GLB_param.th24;
    parameters[6] = GLB_param.th34;
    parameters[7] = GLB_param.th15;
    parameters[8] = GLB_param.th25;
    parameters[9] = GLB_param.th35;
    parameters[10] = GLB_param.th45;
    parameters[11] = GLB_param.th16;
    parameters[12] = GLB_param.th26;
    parameters[13] = GLB_param.th36;
    parameters[14] = GLB_param.th46;
    parameters[15] = GLB_param.th56;
    parameters[16] = GLB_param.dm21sq;
    parameters[17] = GLB_param.dm31sq;
    parameters[18] = GLB_param.dm41sq;
    parameters[19] = GLB_param.dm51sq;
    parameters[20] = GLB_param.dm61sq;
    parameters[21] = GLB_param.delta1;
    parameters[22] = GLB_param.delta2;
    parameters[23] = GLB_param.delta3;

    MLPutRealList(stdlink, parameters, paramnum);
    return;
}

void CalNeuOscGSLML(int ineu,double E,int body_id,double body_params[],long body_params_len,double track_params[],long track_params_len,int neutype,double abs_error,double rel_error,int return_state,int optimization){
    int numneu = GLB_param.numneu;
    double osc_prob[numneu];
    Body body;
    Track track;
    
    if (body_id == 0){
    //body_id = 0 : Vacuum
        // checking parameters lengths
        if (body_params_len == 0){
            Vacuum body;    
        } else {
            MLEvaluate(stdlink, "Message[CalNeuOscGSLML::Invalid body arguments.]");
            MLClearError(stdlink);
            MLNextPacket(stdlink);
            MLNewPacket(stdlink);
            MLPutSymbol(stdlink, "$Failed");
            return;
        }
        
        if (track_params_len == 2){
            Vacuum::Track trackInit(track_params[0],track_params[1]);
            track = trackInit;
        } else {
            MLEvaluate(stdlink, "Message[CalNeuOscGSLML::Invalid track arguments.]");
            MLClearError(stdlink);
            MLNextPacket(stdlink);
            MLNewPacket(stdlink);
            MLPutSymbol(stdlink, "$Failed");
            return;
        }
    } else if(body_id == 1){
    // body_id = 1 : ConstantDensity
        // checking parameters lenghts
        if (body_params_len == 2){
            ConstantDensity body(body_params[0],body_params[1]);    
        } else {
            MLEvaluate(stdlink, "Message[CalNeuOscGSLML::Invalid body arguments.]");
            MLClearError(stdlink);
            MLNextPacket(stdlink);
            MLNewPacket(stdlink);
            MLPutSymbol(stdlink, "$Failed");
            return;
        }
        
        if (track_params_len == 2){
            ConstantDensity::Track trackInit(track_params[0],track_params[1]);
            track = trackInit;
        } else {
            MLEvaluate(stdlink, "Message[CalNeuOscGSLML::Invalid track arguments.]");
            MLClearError(stdlink);
            MLNextPacket(stdlink);
            MLNewPacket(stdlink);
            MLPutSymbol(stdlink, "$Failed");
            return;
        }    
    } else if(body_id == 2){
    // body_id = 2 : VariableDensity
        // checking parameters lenghts
        if (body_params_len == 3){
            VariableDensity body(body_params[0],body_params[1],body_params[2]);    
        } else {
            MLEvaluate(stdlink, "Message[CalNeuOscGSLML::Invalid body arguments.]");
            MLClearError(stdlink);
            MLNextPacket(stdlink);
            MLNewPacket(stdlink);
            MLPutSymbol(stdlink, "$Failed");
            return;
        }
        
        if (track_params_len == 2){
            VariableDensity::Track trackInit(track_params[0],track_params[1]);
            track = trackInit;
        } else {
            MLEvaluate(stdlink, "Message[CalNeuOscGSLML::Invalid track arguments.]");
            MLClearError(stdlink);
            MLNextPacket(stdlink);
            MLNewPacket(stdlink);
            MLPutSymbol(stdlink, "$Failed");
            return;
        }    
    } else if(body_id == 3){
    // body_id = 2 : Star
        // checking parameters lenghts
        if (body_params_len == 3){
            Star body(body_params[0],body_params[1],body_params[2]);    
        } else {
            MLEvaluate(stdlink, "Message[CalNeuOscGSLML::Invalid body arguments.]");
            MLClearError(stdlink);
            MLNextPacket(stdlink);
            MLNewPacket(stdlink);
            MLPutSymbol(stdlink, "$Failed");
            return;
        }
        
        if (track_params_len == 2){
            Star::Track trackInit(track_params[0],track_params[1]);
            track = trackInit;
        } else {
            MLEvaluate(stdlink, "Message[CalNeuOscGSLML::Invalid track arguments.]");
            MLClearError(stdlink);
            MLNextPacket(stdlink);
            MLNewPacket(stdlink);
            MLPutSymbol(stdlink, "$Failed");
            return;
        }    
    } else if(body_id == 4){
    // body_id = 4 : Earth
        // checking parameters lenghts
        if (body_params_len == 0){
            Earth body();    
        } else {
            MLEvaluate(stdlink, "Message[CalNeuOscGSLML::Invalid body arguments.]");
            MLClearError(stdlink);
            MLNextPacket(stdlink);
            MLNewPacket(stdlink);
            MLPutSymbol(stdlink, "$Failed");
            return;
        }
        
        if (track_params_len == 3){
            Earth::Track trackInit(track_params[0],track_params[1],track_params[2]);
            track = trackInit;
        } else {
            MLEvaluate(stdlink, "Message[CalNeuOscGSLML::Invalid track arguments.]");
            MLClearError(stdlink);
            MLNextPacket(stdlink);
            MLNewPacket(stdlink);
            MLPutSymbol(stdlink, "$Failed");
            return;
        }    
    } else if(body_id == 5){
    // body_id = 5 : AGN
        // checking parameters lenghts
        if (body_params_len == 3){
            AGN body(body_params[0],body_params[1],body_params[2]);    
        } else {
            MLEvaluate(stdlink, "Message[CalNeuOscGSLML::Invalid body arguments.]");
            MLClearError(stdlink);
            MLNextPacket(stdlink);
            MLNewPacket(stdlink);
            MLPutSymbol(stdlink, "$Failed");
            return;
        }
        
        if (track_params_len == 2){
            AGN::Track trackInit(track_params[0],track_params[1]);
            track = trackInit;
        } else {
            MLEvaluate(stdlink, "Message[CalNeuOscGSLML::Invalid track arguments.]");
            MLClearError(stdlink);
            MLNextPacket(stdlink);
            MLNewPacket(stdlink);
            MLPutSymbol(stdlink, "$Failed");
            return;
        }    
    } else {
        MLEvaluate(stdlink, "Message[CalNeuOscGSLML::Invalid body name.]");
        MLClearError(stdlink);
        MLNextPacket(stdlink);
        MLNewPacket(stdlink);
        MLPutSymbol(stdlink, "$Failed");
        return;
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
        MLPutRealList(stdlink, osc_prob, 2*numneu);
    } else {
        MLPutRealList(stdlink, osc_prob, numneu);
    }
    
    return;
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
    status = MLMain(argc,argv);
    
    return status;
}


