/*
 * This file automatically produced by /usr/local/Wolfram/Mathematica/8.0/SystemFiles/Links/MathLink/DeveloperKit/Linux-x86-64/CompilerAdditions/mprep from:
 *	neutrinocommon.tm
 * mprep Revision 16 Copyright (c) Wolfram Research, Inc. 1990-2009
 */

#define MPREP_REVISION 16

#include "mathlink.h"

int MLAbort = 0;
int MLDone  = 0;
long MLSpecialCharacter = '\0';

MLINK stdlink = 0;
MLEnvironment stdenv = 0;
#if MLINTERFACE >= 3
MLYieldFunctionObject stdyielder = (MLYieldFunctionObject)0;
MLMessageHandlerObject stdhandler = (MLMessageHandlerObject)0;
#else
MLYieldFunctionObject stdyielder = 0;
MLMessageHandlerObject stdhandler = 0;
#endif /* MLINTERFACE >= 3 */

/********************************* end header *********************************/


# line 1 "neutrinocommon.tm"
// -------------------------------------------
// neutrinocommon.tm : MathLink template file
// -------------------------------------------

// Function declarations

// Physics Parameters Management

# line 38 "neutrinocommon-tm.c"


# line 33 "neutrinocommon.tm"
// Oscillation Calculation

# line 44 "neutrinocommon-tm.c"


# line 44 "neutrinocommon.tm"
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


# line 402 "neutrinocommon-tm.c"


int SetDefaultParameters P(( void));

#if MLPROTOTYPES
static int _tr0( MLINK mlp)
#else
static int _tr0(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	int _tp0;
	if ( ! MLNewPacket(mlp) ) goto L0;

	_tp0 = SetDefaultParameters();

	res = MLAbort ?
		MLPutFunction( mlp, "Abort", 0) : MLPutInteger( mlp, _tp0);

L0:	return res;
} /* _tr0 */


int SetParamsML P(( int _tp1, double _tp2));

#if MLPROTOTYPES
static int _tr1( MLINK mlp)
#else
static int _tr1(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	int _tp0;
	int _tp1;
	double _tp2;
	if ( ! MLGetInteger( mlp, &_tp1) ) goto L0;
	if ( ! MLGetReal( mlp, &_tp2) ) goto L1;
	if ( ! MLNewPacket(mlp) ) goto L2;

	_tp0 = SetParamsML(_tp1, _tp2);

	res = MLAbort ?
		MLPutFunction( mlp, "Abort", 0) : MLPutInteger( mlp, _tp0);
L2: L1: 
L0:	return res;
} /* _tr1 */


void GetParamsML P(( void));

#if MLPROTOTYPES
static int _tr2( MLINK mlp)
#else
static int _tr2(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	if ( ! MLNewPacket(mlp) ) goto L0;
	if( !mlp) return res; /* avoid unused parameter warning */

	GetParamsML();

	res = 1;

L0:	return res;
} /* _tr2 */


void CalNeuOscGSLML P(( int _tp1, double _tp2, int _tp3, double * _tp4, long _tpl4, double * _tp5, long _tpl5, int _tp6, double _tp7, double _tp8, int _tp9, int _tp10));

#if MLPROTOTYPES
static int _tr3( MLINK mlp)
#else
static int _tr3(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	int _tp1;
	double _tp2;
	int _tp3;
	double * _tp4;
	long _tpl4;
	double * _tp5;
	long _tpl5;
	int _tp6;
	double _tp7;
	double _tp8;
	int _tp9;
	int _tp10;
	if ( ! MLGetInteger( mlp, &_tp1) ) goto L0;
	if ( ! MLGetReal( mlp, &_tp2) ) goto L1;
	if ( ! MLGetInteger( mlp, &_tp3) ) goto L2;
	if ( ! MLGetRealList( mlp, &_tp4, &_tpl4) ) goto L3;
	if ( ! MLGetRealList( mlp, &_tp5, &_tpl5) ) goto L4;
	if ( ! MLGetInteger( mlp, &_tp6) ) goto L5;
	if ( ! MLGetReal( mlp, &_tp7) ) goto L6;
	if ( ! MLGetReal( mlp, &_tp8) ) goto L7;
	if ( ! MLGetInteger( mlp, &_tp9) ) goto L8;
	if ( ! MLGetInteger( mlp, &_tp10) ) goto L9;
	if ( ! MLNewPacket(mlp) ) goto L10;

	CalNeuOscGSLML(_tp1, _tp2, _tp3, _tp4, _tpl4, _tp5, _tpl5, _tp6, _tp7, _tp8, _tp9, _tp10);

	res = 1;
L10: L9: L8: L7: L6: L5:	MLDisownRealList( mlp, _tp5, _tpl5);
L4:	MLDisownRealList( mlp, _tp4, _tpl4);
L3: L2: L1: 
L0:	return res;
} /* _tr3 */


static struct func {
	int   f_nargs;
	int   manual;
	int   (*f_func)P((MLINK));
	const char  *f_name;
	} _tramps[4] = {
		{ 0, 0, _tr0, "SetDefaultParameters" },
		{ 2, 0, _tr1, "SetParamsML" },
		{ 0, 0, _tr2, "GetParamsML" },
		{10, 0, _tr3, "CalNeuOscGSLML" }
		};

#define CARDOF_EVALSTRS 0

static int _definepattern P(( MLINK, char*, char*, int));

int  _MLDoCallPacket P(( MLINK, struct func[], int));


#if MLPROTOTYPES
int MLInstall( MLINK mlp)
#else
int MLInstall(mlp) MLINK mlp;
#endif
{
	int _res;
	_res = MLConnect(mlp);
	if (_res) _res = _definepattern(mlp, (char *)"SetDefaultParameters[]", (char *)"{}", 0);
	if (_res) _res = _definepattern(mlp, (char *)"SetParamsML[name_Integer,value_?NumericQ]", (char *)"{name,value}", 1);
	if (_res) _res = _definepattern(mlp, (char *)"GetParamsML[]", (char *)"{}", 2);
	if (_res) _res = _definepattern(mlp, (char *)"CalNeuOscGSLML[ineu_Integer,energy_?NumericQ,bodyname_Integer,bodyparams_List,trackparams_List,neutype_Integer,abserror_?NumericQ,relerror_?NumericQ,returnstate_Integer,optimization_Integer]", (char *)"{ineu,energy,bodyname,bodyparams,trackparams,neutype,abserror,relerror,returnstate,optimization}", 3);
	if (_res) _res = MLPutSymbol( mlp, "End");
	if (_res) _res = MLFlush( mlp);
	return _res;
} /* MLInstall */


#if MLPROTOTYPES
int MLDoCallPacket( MLINK mlp)
#else
int MLDoCallPacket( mlp) MLINK mlp;
#endif
{
	return _MLDoCallPacket( mlp, _tramps, 4);
} /* MLDoCallPacket */

/******************************* begin trailer ********************************/

#ifndef EVALSTRS_AS_BYTESTRINGS
#	define EVALSTRS_AS_BYTESTRINGS 1
#endif


#if CARDOF_EVALSTRS
#if MLPROTOTYPES
static int  _doevalstr( MLINK mlp, int n)
#else
static int  _doevalstr( mlp, n)
	 MLINK mlp; int n;
#endif
{
	long bytesleft, charsleft, bytesnow;
#if !EVALSTRS_AS_BYTESTRINGS
	long charsnow;
#endif
	char **s, **p;
	char *t;

	s = (char **)evalstrs;
	while( n-- > 0){
		if( *s == 0) break;
		while( *s++ != 0){}
	}
	if( *s == 0) return 0;
	bytesleft = 0;
	charsleft = 0;
	p = s;
	while( *p){
		t = *p; while( *t) ++t;
		bytesnow = t - *p;
		bytesleft += bytesnow;
		charsleft += bytesnow;
#if !EVALSTRS_AS_BYTESTRINGS
		t = *p;
		charsleft -= MLCharacterOffset( &t, t + bytesnow, bytesnow);
		/* assert( t == *p + bytesnow); */
#endif
		++p;
	}


	MLPutNext( mlp, MLTKSTR);
#if EVALSTRS_AS_BYTESTRINGS
	p = s;
	while( *p){
		t = *p; while( *t) ++t;
		bytesnow = t - *p;
		bytesleft -= bytesnow;
		MLPut8BitCharacters( mlp, bytesleft, (unsigned char*)*p, bytesnow);
		++p;
	}
#else
	MLPut7BitCount( mlp, charsleft, bytesleft);
	p = s;
	while( *p){
		t = *p; while( *t) ++t;
		bytesnow = t - *p;
		bytesleft -= bytesnow;
		t = *p;
		charsnow = bytesnow - MLCharacterOffset( &t, t + bytesnow, bytesnow);
		/* assert( t == *p + bytesnow); */
		charsleft -= charsnow;
		MLPut7BitCharacters(  mlp, charsleft, *p, bytesnow, charsnow);
		++p;
	}
#endif
	return MLError( mlp) == MLEOK;
}
#endif /* CARDOF_EVALSTRS */


#if MLPROTOTYPES
static int  _definepattern( MLINK mlp, char *patt, char *args, int func_n)
#else
static int  _definepattern( mlp, patt, args, func_n)
	MLINK  mlp;
	char  *patt, *args;
	int    func_n;
#endif
{
	MLPutFunction( mlp, "DefineExternal", (long)3);
	  MLPutString( mlp, patt);
	  MLPutString( mlp, args);
	  MLPutInteger( mlp, func_n);
	return !MLError(mlp);
} /* _definepattern */


#if MLPROTOTYPES
int _MLDoCallPacket( MLINK mlp, struct func functable[], int nfuncs)
#else
int _MLDoCallPacket( mlp, functable, nfuncs)
	MLINK mlp;
	struct func functable[];
	int nfuncs;
#endif
{
	long len;
	int n, res = 0;
	struct func* funcp;

	if( ! MLGetInteger( mlp, &n) ||  n < 0 ||  n >= nfuncs) goto L0;
	funcp = &functable[n];

	if( funcp->f_nargs >= 0
	&& ( ! MLCheckFunction(mlp, "List", &len)
	     || ( !funcp->manual && (len != funcp->f_nargs))
	     || (  funcp->manual && (len <  funcp->f_nargs))
	   )
	) goto L0;

	stdlink = mlp;
	res = (*funcp->f_func)( mlp);

L0:	if( res == 0)
		res = MLClearError( mlp) && MLPutSymbol( mlp, "$Failed");
	return res && MLEndPacket( mlp) && MLNewPacket( mlp);
} /* _MLDoCallPacket */


#if MLPROTOTYPES
mlapi_packet MLAnswer( MLINK mlp)
#else
mlapi_packet MLAnswer( mlp)
	MLINK mlp;
#endif
{
	mlapi_packet pkt = 0;

	while( !MLDone && !MLError(mlp) && (pkt = MLNextPacket(mlp), pkt) && pkt == CALLPKT){
		MLAbort = 0;
		if( !MLDoCallPacket(mlp)) pkt = 0;
	}
	MLAbort = 0;
	return pkt;
} /* MLAnswer */



/*
	Module[ { me = $ParentLink},
		$ParentLink = contents of RESUMEPKT;
		Message[ MessageName[$ParentLink, "notfe"], me];
		me]
*/

#if MLPROTOTYPES
static int refuse_to_be_a_frontend( MLINK mlp)
#else
static int refuse_to_be_a_frontend( mlp)
	MLINK mlp;
#endif
{
	int pkt;

	MLPutFunction( mlp, "EvaluatePacket", 1);
	  MLPutFunction( mlp, "Module", 2);
	    MLPutFunction( mlp, "List", 1);
		  MLPutFunction( mlp, "Set", 2);
		    MLPutSymbol( mlp, "me");
	        MLPutSymbol( mlp, "$ParentLink");
	  MLPutFunction( mlp, "CompoundExpression", 3);
	    MLPutFunction( mlp, "Set", 2);
	      MLPutSymbol( mlp, "$ParentLink");
	      MLTransferExpression( mlp, mlp);
	    MLPutFunction( mlp, "Message", 2);
	      MLPutFunction( mlp, "MessageName", 2);
	        MLPutSymbol( mlp, "$ParentLink");
	        MLPutString( mlp, "notfe");
	      MLPutSymbol( mlp, "me");
	    MLPutSymbol( mlp, "me");
	MLEndPacket( mlp);

	while( (pkt = MLNextPacket( mlp), pkt) && pkt != SUSPENDPKT)
		MLNewPacket( mlp);
	MLNewPacket( mlp);
	return MLError( mlp) == MLEOK;
}


#if MLPROTOTYPES
#if MLINTERFACE >= 3
int MLEvaluate( MLINK mlp, char *s)
#else
int MLEvaluate( MLINK mlp, charp_ct s)
#endif /* MLINTERFACE >= 3 */
#else
int MLEvaluate( mlp, s)
	MLINK mlp;
#if MLINTERFACE >= 3
	char *s;
#else
	charp_ct s;
#endif /* MLINTERFACE >= 3 */
#endif
{
	if( MLAbort) return 0;
	return MLPutFunction( mlp, "EvaluatePacket", 1L)
		&& MLPutFunction( mlp, "ToExpression", 1L)
		&& MLPutString( mlp, s)
		&& MLEndPacket( mlp);
} /* MLEvaluate */


#if MLPROTOTYPES
#if MLINTERFACE >= 3
int MLEvaluateString( MLINK mlp, char *s)
#else
int MLEvaluateString( MLINK mlp, charp_ct s)
#endif /* MLINTERFACE >= 3 */
#else
int MLEvaluateString( mlp, s)
	MLINK mlp;
#if MLINTERFACE >= 3
	char *s;
#else
	charp_ct s;
#endif /* MLINTERFACE >= 3 */
#endif
{
	int pkt;
	if( MLAbort) return 0;
	if( MLEvaluate( mlp, s)){
		while( (pkt = MLAnswer( mlp), pkt) && pkt != RETURNPKT)
			MLNewPacket( mlp);
		MLNewPacket( mlp);
	}
	return MLError( mlp) == MLEOK;
} /* MLEvaluateString */


#if MLINTERFACE >= 3
#if MLPROTOTYPES
void MLDefaultHandler( MLINK mlp, int message, int n)
#else
void MLDefaultHandler( mlp, message, n)
	MLINK mlp;
	int message, n;
#endif
#else
#if MLPROTOTYPES
void MLDefaultHandler( MLINK mlp, unsigned long message, unsigned long n)
#else
void MLDefaultHandler( mlp, message, n)
	MLINK mlp;
	unsigned long message, n;
#endif
#endif /* MLINTERFACE >= 3 */
{
	switch (message){
	case MLTerminateMessage:
		MLDone = 1;
	case MLInterruptMessage:
	case MLAbortMessage:
		MLAbort = 1;
	default:
		return;
	}
}

#if MLPROTOTYPES
#if MLINTERFACE >= 3
static int _MLMain( char **argv, char **argv_end, char *commandline)
#else
static int _MLMain( charpp_ct argv, charpp_ct argv_end, charp_ct commandline)
#endif /* MLINTERFACE >= 3 */
#else
static int _MLMain( argv, argv_end, commandline)
#if MLINTERFACE >= 3
  char **argv, argv_end;
  char *commandline;
#else
  charpp_ct argv, argv_end;
  charp_ct commandline;
#endif /* MLINTERFACE >= 3 */
#endif
{
	MLINK mlp;
#if MLINTERFACE >= 3
	int err;
#else
	long err;
#endif /* MLINTERFACE >= 3 */

	if( !stdenv)
		stdenv = MLInitialize( (MLParametersPointer)0);
	if( stdenv == (MLEnvironment)0) goto R0;

#if MLINTERFACE >= 3
	if( !stdhandler)
		stdhandler = (MLMessageHandlerObject)MLDefaultHandler;
#else
	if( !stdhandler)
		stdhandler = MLCreateMessageHandler( stdenv, MLDefaultHandler, 0);
#endif /* MLINTERFACE >= 3 */


	mlp = commandline
		? MLOpenString( stdenv, commandline, &err)
#if MLINTERFACE >= 3
		: MLOpenArgcArgv( stdenv, (int)(argv_end - argv), argv, &err);
#else
		: MLOpenArgv( stdenv, argv, argv_end, &err);
#endif
	if( mlp == (MLINK)0){
		MLAlert( stdenv, MLErrorString( stdenv, err));
		goto R1;
	}

	if( stdyielder) MLSetYieldFunction( mlp, stdyielder);
	if( stdhandler) MLSetMessageHandler( mlp, stdhandler);

	if( MLInstall( mlp))
		while( MLAnswer( mlp) == RESUMEPKT){
			if( ! refuse_to_be_a_frontend( mlp)) break;
		}

	MLClose( mlp);
R1:	MLDeinitialize( stdenv);
	stdenv = (MLEnvironment)0;
R0:	return !MLDone;
} /* _MLMain */


#if MLPROTOTYPES
#if MLINTERFACE >= 3
int MLMainString( char *commandline)
#else
int MLMainString( charp_ct commandline)
#endif /* MLINTERFACE >= 3 */
#else
#if MLINTERFACE >= 3
int MLMainString( commandline)  char *commandline;
#else
int MLMainString( commandline)  charp_ct commandline;
#endif /* MLINTERFACE >= 3 */
#endif
{
	return _MLMain( (charpp_ct)0, (charpp_ct)0, commandline);
}

#if MLPROTOTYPES
int MLMainArgv( char** argv, char** argv_end) /* note not FAR pointers */
#else
int MLMainArgv( argv, argv_end)  char **argv, **argv_end;
#endif
{   
	static char FAR * far_argv[128];
	int count = 0;
	
	while(argv < argv_end)
		far_argv[count++] = *argv++;
		 
	return _MLMain( far_argv, far_argv + count, (charp_ct)0);

}

#if MLPROTOTYPES
#if MLINTERFACE >= 3
int MLMain( int argc, char **argv)
#else
int MLMain( int argc, charpp_ct argv)
#endif /* MLINTERFACE >= 3 */
#else
#if MLINTERFACE >= 3
int MLMain( argc, argv) int argc; char **argv;
#else
int MLMain( argc, argv) int argc; charpp_ct argv;
#endif /* MLINTERFACE >= 3 */
#endif
{
#if MLINTERFACE >= 3
 	return _MLMain( argv, argv + argc, (char *)0);
#else
 	return _MLMain( argv, argv + argc, (charp_ct)0);
#endif /* MLINTERFACE >= 3 */
}
 
