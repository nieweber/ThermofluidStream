within ThermofluidStream.Processes.Internal.FlowResistance;
function zetaPressureLoss "Pressure loss function based on zeta value"
  extends Internal.FlowResistance.partialPressureLoss;

  input ThermofluidStream.Processes.Internal.GeometryOfResistance geometry = ThermofluidStream.Processes.Internal.GeometryOfResistance.circular
  "Geometry of cross sectional area"
    annotation(Dialog(enable=true),
     choices(
      choice=ThermofluidStream.Processes.Internal.GeometryOfResistance.circular "Circular",
      choice=ThermofluidStream.Processes.Internal.GeometryOfResistance.rectangle "Rectangle",
      choice=ThermofluidStream.Processes.Internal.GeometryOfResistance.other "Other"));

  input SI.Length r(min=0) "Pipe radius" annotation(Dialog(enable = (geometry == ThermofluidStream.Processes.Internal.GeometryOfResistance.circular)));
  input SI.Length a(min=0) = 0 "Rectangle width"
    annotation(Dialog(enable = (geometry == ThermofluidStream.Processes.Internal.GeometryOfResistance.rectangle)));
  input SI.Length b(min=0) = 0 "Rectangle height"
    annotation(Dialog(enable = (geometry == ThermofluidStream.Processes.Internal.GeometryOfResistance.rectangle)));

  input SI.Area A = Modelica.Constants.pi*r*r "Reference area from parameter"
    annotation(Dialog(enable=(geometry == ThermofluidStream.Processes.Internal.GeometryOfResistance.other)));

  input Real zeta( unit = "1") "Zeta value of component"
    annotation(Dialog(enable=true));

  output Real zeta_value;

protected
  SI.Area A_zeta "Reference area either from radius or set by parameter";

algorithm

  if geometry == ThermofluidStream.Processes.Internal.GeometryOfResistance.circular then
    A_zeta :=Modelica.Constants.pi*r*r;
  elseif geometry == ThermofluidStream.Processes.Internal.GeometryOfResistance.rectangle then
    A_zeta := a*b;
  elseif geometry == ThermofluidStream.Processes.Internal.GeometryOfResistance.other then
    A_zeta := A;
  end if;

  pressureLoss :=zeta/(2*rho)*Modelica.Fluid.Utilities.regSquare(m_flow/A_zeta);

  zeta_value:=zeta;

  annotation (Documentation(info="<html>
<p>For specific components (armatures, fittings, pipe sections, grids, ...), the zeta value is often given in the data sheet.</p>
<p>Together with a given reference area A, the pressure drop can be calculated:</p>
<p style=\"margin-left: 40px;\"><span style=\"font-family: Courier New;\">dp := zeta/(2*rho)*m_flow^2/A^2</span></p>
<p>The square of the mass-flow is regularized using the regSquare function from the MSL.</p>
</html>"));
end zetaPressureLoss;
