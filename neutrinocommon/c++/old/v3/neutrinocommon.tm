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
:Function:      SetParams
:Pattern:       SetParams[paramname_String,value_?NumericQ]
:Arguments:     {name,value}
:ArgumentTypes: {String,Real}
:ReturnType:    Integer
:End:

:Begin:
:Function:      GetParams
:Pattern:       GetParams[]
:Arguments:     {}
:ArgumentTypes: {}
:ReturnType:    Manual
:End:

// Oscillation Calculation

:Begin:
:Function:      CalNeuOscGSLML
:Pattern:       CalNeuOscGSLML[ineu_Integer,energy_Real,bodyname_String,bodyparams_RealList,trackparams_RealList,neutype_String,abserror_Real,relerror_Real,returnstate_Integer,optimization_Integer]
:Arguments:     {ineu,E,bodyname,bodyparams,track_params,neutype,abserror,relerror,returnstate,optimization}
:ArgumentTypes: {Integer,Real,String,RealList,RealList,String,Real,Real,Integer,Integer}
:ReturnType:    Manual
:End:


