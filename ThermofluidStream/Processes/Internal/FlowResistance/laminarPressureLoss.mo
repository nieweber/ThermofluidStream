within ThermofluidStream.Processes.Internal.FlowResistance;
function laminarPressureLoss
  "Laminar flow loss function (Hagen-Poiseuille)"
  extends Internal.FlowResistance.partialPressureLoss;

  input ShapeOfResistance shape=ThermofluidStream.Processes.Internal.ShapeOfResistance.circular "Shape of cross sectional area"
    annotation (Dialog(enable=true),
    choices(
      choice=ThermofluidStream.Processes.Internal.ShapeOfResistance.circular "Circular",
      choice=ThermofluidStream.Processes.Internal.ShapeOfResistance.rectangle "Rectangle",
      choice=ThermofluidStream.Processes.Internal.ShapeOfResistance.other "Other"));

  input SI.Length l(min=0) "Length of component" annotation(Dialog(enable=true));

  input SI.Length r(min=0) "Pipe radius" annotation(Dialog(enable=(shape == ThermofluidStream.Processes.Internal.ShapeOfResistance.circular)));

  input SI.Length a(min=0) = 0 "Rectangle width"
    annotation(Dialog(enable=(shape == ThermofluidStream.Processes.Internal.ShapeOfResistance.rectangle)));
  input SI.Length b(min=0) = 0 "Rectangle height"
    annotation(Dialog(enable=(shape == ThermofluidStream.Processes.Internal.ShapeOfResistance.rectangle)));
  input SI.Length d_h_input = 0 "Custom hydraulic diameter if shape not available" annotation (Dialog(enable=(shape == ThermofluidStream.Processes.Internal.ShapeOfResistance.other)));

  import Modelica.Constants.pi;

algorithm

  if shape == ThermofluidStream.Processes.Internal.ShapeOfResistance.circular then
    d_h := 2*r;
  elseif shape == ThermofluidStream.Processes.Internal.ShapeOfResistance.rectangle then
    d_h := 2*a*b/(a+b);
  elseif shape == ThermofluidStream.Processes.Internal.ShapeOfResistance.other then
    d_h := d_h_input;
  end if;

  if not shape == ThermofluidStream.Processes.Internal.ShapeOfResistance.circular then
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
