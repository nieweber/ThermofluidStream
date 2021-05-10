within ThermofluidStream.Boundaries.Internal;
type InitializationMethodsPhaseSeperator = enumeration(
    h "specific enthalpy h_0",
    M "Mass M_0",
    l "Liquid Level l_0",
    x "Vapor Quality x_0") "Choices for Initialization of state h of PhaseSeperator"
  annotation (Icon(coordinateSystem(preserveAspectRatio=false)),
    Diagram(coordinateSystem(preserveAspectRatio=false)),
    Documentation(info="<html>
        <p><span style=\"font-family: Courier New;\">Choices for Initailtaion of a state h.</span></p>
        </html>"));
