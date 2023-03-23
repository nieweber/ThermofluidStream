within ThermofluidStream.Processes.Internal.FlowResistance;
function laminarPressureLoss
  "Laminar flow loss function (Hagen-Poiseuille)"
  extends Internal.FlowResistance.partialPressureLoss;

  input ThermofluidStream.Processes.Internal.GeometryOfResistance geometry = ThermofluidStream.Processes.Internal.GeometryOfResistance.circular
  "Geometry of cross sectional area"
    annotation(Dialog(enable=true),
     choices(
      choice=ThermofluidStream.Processes.Internal.GeometryOfResistance.circular "Circular",
      choice=ThermofluidStream.Processes.Internal.GeometryOfResistance.rectangle "Rectangle",
      choice=ThermofluidStream.Processes.Internal.GeometryOfResistance.other "Other"));

  input SI.Length l(min=0) "Length of component" annotation(Dialog(enable=true));

  input SI.Length r(min=0) "Pipe radius" annotation(Dialog(enable = (geometry == ThermofluidStream.Processes.Internal.GeometryOfResistance.circular)));

  input SI.Length a(min=0) = 0 "Rectangle width"
    annotation(Dialog(enable = (geometry == ThermofluidStream.Processes.Internal.GeometryOfResistance.rectangle)));
  input SI.Length b(min=0) = 0 "Rectangle height"
    annotation(Dialog(enable = (geometry == ThermofluidStream.Processes.Internal.GeometryOfResistance.rectangle)));
  input SI.Length d_h_input = 0 "Custom hydraulic diameter if shape not available" annotation (Dialog(enable = (geometry == ThermofluidStream.Processes.Internal.GeometryOfResistance.other)));

  import Modelica.Constants.pi;

algorithm

  if geometry == ThermofluidStream.Processes.Internal.GeometryOfResistance.circular then
    d_h := 2*r;
  elseif geometry == ThermofluidStream.Processes.Internal.GeometryOfResistance.rectangle then
    d_h := 2*a*b/(a+b);
  elseif geometry == ThermofluidStream.Processes.Internal.GeometryOfResistance.other then
    d_h := d_h_input;
  end if;

  if  not geometry == ThermofluidStream.Processes.Internal.GeometryOfResistance.circular then
    r_h :=d_h/4;
  else
    r_h :=r;
  end if;

  pressureLoss := m_flow * (8*mu*l)/(pi*rho*r_h^4);

  annotation (Documentation(info="<html>
<p>
Pressure loss after Hagen-Poiseuille:
</p>
<blockquote><pre>
pressureLoss := m_flow * (8*mu*l)/(pi*rho*r^4);
</pre></blockquote>
</html>"));
end laminarPressureLoss;
