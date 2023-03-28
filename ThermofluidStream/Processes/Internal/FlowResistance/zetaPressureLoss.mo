within ThermofluidStream.Processes.Internal.FlowResistance;
function zetaPressureLoss "Pressure loss function based on zeta value"
  extends Internal.FlowResistance.partialPressureLoss;

  input ShapeOfResistance shape=ThermofluidStream.Processes.Internal.ShapeOfResistance.circular "Shape of cross sectional area"
  annotation (Dialog(enable=true),
  choices(
      choice=ThermofluidStream.Processes.Internal.ShapeOfResistance.circular "Circular",
      choice=ThermofluidStream.Processes.Internal.ShapeOfResistance.rectangle "Rectangle",
      choice=ThermofluidStream.Processes.Internal.ShapeOfResistance.other "Other"));

  input SI.Length r(min=0) "Pipe radius" annotation(Dialog(enable=(shape == ThermofluidStream.Processes.Internal.ShapeOfResistance.circular)));
  input SI.Length a(min=0) = 0 "Rectangle width"
    annotation(Dialog(enable=(shape == ThermofluidStream.Processes.Internal.ShapeOfResistance.rectangle)));
  input SI.Length b(min=0) = 0 "Rectangle height"
    annotation(Dialog(enable=(shape == ThermofluidStream.Processes.Internal.ShapeOfResistance.rectangle)));

  input SI.Area A = Modelica.Constants.pi*r*r "Reference area from parameter"
    annotation(Dialog(enable=(shape == ThermofluidStream.Processes.Internal.ShapeOfResistance.other)));

  input Real zeta( unit = "1") "Zeta value of component"
    annotation(Dialog(enable=true));

  input Boolean compressible = false "Use pressure loss function for compressible Media?"
    annotation (Dialog(enable = true));


protected
  SI.Area A_zeta "Reference area either from radius or set by parameter";

algorithm

  if shape == ThermofluidStream.Processes.Internal.ShapeOfResistance.circular then
    d_h := 2*r;
    A_zeta :=Modelica.Constants.pi*r*r;
  elseif shape == ThermofluidStream.Processes.Internal.ShapeOfResistance.rectangle then
    d_h := 2*a*b/(a+b);
    A_zeta := Modelica.Constants.pi*d_h/4;
  elseif shape == ThermofluidStream.Processes.Internal.ShapeOfResistance.other then
    A_zeta := A;
  end if;

  if compressible then
    result.dp := (p_out - (p_out*p_out+zeta*Rs*T*Modelica.Fluid.Utilities.regSquare(m_flow/A_zeta))^0.5)*(-1);
  else
    result.dp := zeta/(2*rho)*Modelica.Fluid.Utilities.regSquare(m_flow/A_zeta);
  end if;

  result.zeta :=zeta;
  result.A := A_zeta;
  result.v := m_flow/(rho*A_zeta);

  annotation (Documentation(info="<html>
<p>For specific components (armatures, fittings, pipe sections, grids, ...), the zeta value is often given in the data sheet.</p>
<p>Together with a given reference area A, the pressure drop can be calculated:</p>
<p style=\"margin-left: 40px;\"><span style=\"font-family: Courier New;\">dp := zeta/(2*rho)*m_flow^2/A^2</span></p>
<p>The square of the mass-flow is regularized using the regSquare function from the MSL.</p>
</html>"));
end zetaPressureLoss;
