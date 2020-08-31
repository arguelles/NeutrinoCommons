// -------------------------------------------
// neutrinocommon.tm : MathLink template file
// -------------------------------------------

// Function declarations

// Physics Parameters Management

:Begin:
:Function:      SetDefaultParameters
:Pattern:       SetDefaultParameters[]
:Arguments:     {}
:ArgumentsTypes:{}
:ReturnType:    Integer
:End:

:Begin:
:Function:      SetParams
:Pattern:       SetParams[paramname_String,value_?NumericQ]
:Arguments:     {name,value}
:ArgumentsTypes:{String,Real}
:ReturnType:    Integer
:End:

:Begin:
:Function:      GetParams
:Pattern:       GetParams[]
:Arguments:     {}
:ArgumentsTypes:{}
:ReturnType:    Manual
:End:

// Oscillation Calculation

:Begin:
:Function:      CalNeuOscGSLML
:Pattern:       CalNeuOscGSLML[ineu_Integer,energy_Real,bodyname_String,bodyparams_RealList,trackparams_RealList,neutype_String,abserror_Real,relerror_Real,returnstate_Boolean,optimization_Boolean]
:Arguments:     {ineu,E,bodyname,bodyparams,track_params,neutype,abserror,relerror,returnstate,optimization}
:ArgumentsTypes:{Integer,Real,String,RealList,RealList,String,Real,Real,Boolean,Boolean}
:ReturnType:    Manual
:End:


