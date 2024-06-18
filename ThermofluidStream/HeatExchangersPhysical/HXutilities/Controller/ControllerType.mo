within ThermofluidStream.HeatExchangersPhysical.HXutilities.Controller;
type ControllerType = enumeration(
    P
     "P controller",
    I
     "I controller",
    PI
      "PI controller",
    PD
      "PD controller",
    PID
       "PID controller") "Enumeration defining P, I, PI, PD, or PID simple controller type"
                                                                 annotation (
  Evaluate=true);
