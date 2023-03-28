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

  input Boolean useMeanDensity = true "Use mean density for dp calculation?" annotation(Evaluate = true, Dialog(enable = true));

  import Modelica.Constants.pi;
protected
  SI.Density rho_mean "Mean density";
  SI.Area hydraulicArea "Hydraulic cross section area";

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

  rho_mean :=0.5*(rho + rho_out);
  hydraulicArea := pi*r_h^2;

  if not useMeanDensity then
  result.dp := m_flow * (8*mu*l)/(pi*rho*r_h^4);
  else
  result.dp := m_flow * (8*mu*l)/(pi*rho_mean*r_h^4);
  //result.dp := m_flow * (8*mu*l)/(pi*rho*r_h^4);
  //result.dp := p_out - sqrt(p_out*p_out+16*mu*Rs*T*m_flow*l/(pi*r_h^4));
  end if;

  result.A :=hydraulicArea;
  result.d_h :=d_h;
  result.v := m_flow/(rho*hydraulicArea);




  annotation (Documentation(info="<html>
<p>
Pressure loss after Hagen-Poiseuille:
</p>
<blockquote><pre>
pressureLoss := m_flow * (8*mu*l)/(pi*rho*r^4);
</pre></blockquote>
</html>"));
end laminarPressureLoss;
