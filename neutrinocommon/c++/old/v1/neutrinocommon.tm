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
:Pattern:       SetParam[paramname_String,value_?NumericQ]
:Arguments:     {name,value}
:ArgumentsTypes:{String,Real}
:ReturnType:    Integer
:End:

:Begin:
:Function:      GetParams
:Pattern:       GetParams[params_List]
:Arguments:     {params}
:ArgumentsTypes:{List}
:ReturnType:    Manual
:End:

:Begin:
:Function:      GetParameters
:Pattern:       GetParameters[params_List]
:Arguments:     {params}
:ArgumentsTypes:{List}
:ReturnType:    Manual
:End:

// Oscillation Calculation

:Begin:
:Function:      CalNeuOscGSLML
:Pattern:       CalNeuOscGSL[ineu_Integer,energy_Real,track_List,body_List,neutype_String]
:Arguments:     {ineu,energy,track,body,neutype}
:ArgumentsTypes:{Integer,Real,List,List,String}
:ReturnType:    Real
:End:


