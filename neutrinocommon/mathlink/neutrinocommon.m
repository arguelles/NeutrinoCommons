(* ::Package:: *)

(* ----------------------- FILE: neutrinocommon.m ------------------------*)
(*                                                                         *)
(* NeutrinoCommon Mathematica Interface                                   *)
(* Author: Carlos Arguelles (carguellesdel@gmail.com)                    *)
(*                                                                         *)
(* ----------------------------------------------------------------------- *)


BeginPackage["NeutrinoCommon`"];

Clear[Vacuum,VacuumTrack,SetParams,GetParams,Parameters,NumberParameter,ParametersRules,InverseParametersRules];

Uninstall[NeutrinoCommonLink];
NeutrinoCommonLink = Install[NotebookDirectory[]<>"mathlink_src/neutrinocommon.exe"];

(***************************** units and constants ********************************)
(* physical constants *)
    GF = 1.16639 10.0^-23;	          (* [eV^-2] Fermi Constant  *)
    Na = 6.0221415 10.0^+23;		    (* [mol cm^-3] Avogadro Number *)
    G  = 6.67300 10.0^-11;              (* [m^3 kg^-1 s^-2] *)
    alpha = 1.0/137.0;                  (* [dimensionless] fine-structure constant  *)
(*unit conversion factors *)
	(*energy*)
    TeV = 1.0 10.0^12;                   (* [eV/TeV]*)
    GeV = 1.0 10.0^9;                    (* [eV/GeV]*)
    MeV = 1.0 10.0^6;                    (* [eV/MeV]*)
    keV = 1.0 10.0^3;                    (* [eV/keV]*)
    Joule = 1/1.60225 10.0^-19;          (* [eV/J]*)
    (* Mass *)
    kg = 5.62 10.0^35;                   (* [eV/kg] *)
    gr = 1.0 10.0^-3*kg;                 (* [eV/g] *)
    (* Time *)
    sec = 1.523 10.0^15;                 (* [eV^-1/s] *)
    hour = 3600.0 sec;                   (* [eV^-1/h] *)
    day = 24.0 hour;                     (* [eV^-1/d]*)
    year = 365.0 day;                    (* [eV^-1/yr]*)
    (* Distance *)
    meter = 5.076 10.0^6;                (* [eV^-1/m]*)
    cm = 1.0 10.0^-2 meter;              (* [eV^-1/cm]*)
    km = 1.0 10.0^3 meter;               (* [eV^-1/km]*)
    fermi = 1.0 10.0^-15 meter;          (* [eV^-1/fm]*)
    angstrom = 1.0 10.0^-10 meter;       (* [eV^-1/A]*)
    AU = 149.60 10.0^9 meter;            (* [eV^-1/AU]*)
    parsec = 3.08568025 10.0^16 meter;   (* [eV^-1/parsec]*)
    (* luminocity *)
    picobarn = 1.0 10.0^-36 cm^2;       (* [eV^-2/pb]*)
    femtobarn = 1.0 10.0^-39 cm^2;      (* [eV^-2/fb]*)
    (* Presure *)
    Pascal = Joule/meter^3;              (* [eV^4/Pa]*)
    hPascal = 100.0*Pascal;               (* [eV^4/hPa]*)
    atm = 101325.0*Pascal;                (* [eV^4/atm]*)
    psi = 6893.0*Pascal;                 (* [eV^4/psi]*)
    (* Temperature *)
    Kelvin = 1/1.1604505 10.0^4;         (* [eV/K]*)
    (* Angle *)
    degree = pi/180.0;                   (* [rad/degree]*)

(********************* WRAPPING AUXILIARY FUNCTIONS ********************************)
CalNeuOsc[istate_,neuE_,{bdname_,bdparams___},{trkparams___},ntype_,abserr_:10.0^-6,relerr_:10.0^-6,rstate_:0,opt_:0]:=NeutrinoCommon`CalNeuOscGSLML[istate,neuE,bdname,{bdparams},{trkparams},ntype,abserr,relerr,rstate,opt];
(* bodies definitions *)
Vacuum[] := {0};
VacuumTrack[xini_,xend_]:={xini,xend};
ConstantDensity[density_,ye_]:={1,density,ye};
ConstantDensityTrack[xini_,xend_]:={xini,xend};
VariableDensity[mindensity_,maxdensity_,ye_]:={2,mindensity,maxdensity,ye};
VariableDensityTrack[xini_,xend_]:={xini,xend};
Star[radius_,iniradius_,inidensity_]:={3,radius,iniradius,inidensity};
StarTrack[xini_,xend_]:={xini,xend};
Earth[]:={4};
EarthTrack[xini_,xend_,baseline_]:={xini,xend,baseline};
AGN[radius_,inidensity_,ye_]:={5,radius,inidensity,ye};
AGNTrack[xini_,xend_]:={xini,xend};
(* oscillation parameters management *)
Parameters = {"NUMNEU","TH12","TH13","TH23","TH14","TH24","TH34","TH15","TH25","TH35","TH45","TH16","TH26","TH36","TH46","TH56","DM21SQ","DM31SQ","DM41SQ","DM51SQ","DM61SQ","DELTA1","DELTA2","DELTA3"};
NumberParameter = Length[Parameters];
ParametersRules = Table[Parameters[[k]]->(k-1),{k,NumberParameter}];
InverseParametersRules = Table[(k-1)->Parameters[[k]],{k,NumberParameter}];

SetParams[paramname_String,value_Real]:=NeutrinoCommon`SetParamsML[paramname/.ParametersRules,value];
SetParams[x__List]:=SetParams[#1[[1]],#1[[2]]]&/@{x};

GetParams[]:=MapIndexed[{(#2[[1]]-1),#1}&,NeutrinoCommon`GetParamsML[]]//.InverseParametersRules;

Begin["`Private`"];
End[];

EndPackage[];
