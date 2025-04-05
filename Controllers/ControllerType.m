classdef ControllerType < Simulink.IntEnumType
    enumeration
        No_Control(0),
        CasADi(1),
        CasADi_freetf(2),
        CasADi_delay(3),
        CasADi_freetf_delay(4),
        CasADi_freetf_freefinal_6DOF(5)
    end
end
