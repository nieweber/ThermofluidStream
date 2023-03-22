within ThermofluidStream.Processes.Internal.FlowResistance;
function laminarPressureLoss
  "Laminar flow loss function (Hagen-Poiseuille)"
  extends Internal.FlowResistance.partialPressureLoss;

  input SI.Length r(min=0) "Pipe radius" annotation(Dialog(enable = true));
  input SI.Length l(min=0) "Pipe length" annotation(Dialog(enable = true));

  import Modelica.Constants.pi;

algorithm
  pressureLoss := m_flow * (8*mu*l)/(pi*rho*r^4);

  annotation (Documentation(info="<html>
<p>
Pressure loss after Hagen-Poiseuille:
</p>
<blockquote><pre>
pressureLoss := m_flow * (8*mu*l)/(pi*rho*r^4);
</pre></blockquote>
</html>"));
end laminarPressureLoss;
