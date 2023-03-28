within ThermofluidStream.Processes.Internal.FlowResistance;
function laminarTurbulentPressureLoss
  "Laminar and turbulent flow regimes pressure loss function (Cheng 2008)"
  extends Internal.FlowResistance.partialPressureLoss;
  import Modelica.Constants.pi;

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

  input SI.Length ks_input(min=1e-7) = 1e-7 "Pipe roughness"
    annotation(Dialog(enable=(material == ThermofluidStream.Processes.Internal.Material.other)));

  input ThermofluidStream.Processes.Internal.Material material=ThermofluidStream.Processes.Internal.Material.other "Material of pipe"
    annotation (Dialog(enable=true),
     choices(
      choice=ThermofluidStream.Processes.Internal.Material.concrete "Concrete ks=5mm",
      choice=ThermofluidStream.Processes.Internal.Material.wood "Wood ks=0.5mm",
      choice=ThermofluidStream.Processes.Internal.Material.castIron "Cast Iron ks=0.25mm",
      choice=ThermofluidStream.Processes.Internal.Material.galvanizedIron "Galvanized Iron ks=0.15mm",
      choice=ThermofluidStream.Processes.Internal.Material.steel "Steel ks=0.059mm",
      choice=ThermofluidStream.Processes.Internal.Material.drawnPipe "Drawn Pipe ks=0.0015mm"));

protected
  constant Real R_laminar_DarcyWeisbach_min(unit="1") = 500 "Minimal Reynolds number to use the general equation. Laminar flow before";
  SI.Length ks "Pipe roughness";

  Real a_factor(unit="1") "Laminar flow factor for the DarcyWeisbach equation (1=laminar flow; 0=turbulent flow)";
  Real b_factor(unit="1") "Turbulent flow factor for DarcyWeisbach equation (1=fully smooth turbulent flow; 0= fully rough turbulent flow)";
  Real lambda_aux(unit="1") "darcy friction factor for DarcyWeisbach equation";

  SI.Velocity u "Median flow velocity";
  Real Re(unit="1") "Reynolds number for flow through the pipe";
  constant Real eps(unit="1") = 0.001;
algorithm
  if material == ThermofluidStream.Processes.Internal.Material.concrete then
    ks :=5e-3;
  elseif material == ThermofluidStream.Processes.Internal.Material.wood then
    ks :=0.5e-3;
  elseif material == ThermofluidStream.Processes.Internal.Material.castIron then
    ks :=0.25e-3;
  elseif material == ThermofluidStream.Processes.Internal.Material.galvanizedIron then
    ks :=0.15e-3;
  elseif material == ThermofluidStream.Processes.Internal.Material.steel then
    ks :=0.059e-3;
  elseif material == ThermofluidStream.Processes.Internal.Material.drawnPipe then
    ks :=0.0015e-3;
  else
    ks :=ks_input;
  end if;

  assert(ks <r_h, "ks must be smaller than r");

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

  u :=m_flow/(r_h^2*pi*rho);
  // add eps to Re to avoid 0^0 error in computation of lambda_aux for Re=0 (-> a=1).
  Re :=abs(u)*rho*2*r_h/mu + eps;

  // cheng 2008. Formulas for Friction Factor in Transitional Regimes. Journal of Hydraulic Engineering.
  a_factor := 1/(1 + (Re/2720)^9);
  b_factor := 1/(1 + (Re/(160*2*r_h/ks))^2);
  //compute lambda_aux = Re*lambda to avoid devision by zero at Re=0 and to avoid if-else
  lambda_aux :=64^a_factor*Re^(1 - a_factor)*((1.8*log10(Re/6.8))^(2*(a_factor - 1)*b_factor)*(2*log10(3.7*2*r_h/ks))^(2*(a_factor - 1)*(1 - b_factor)));

  result.dp := lambda_aux*l*mu*u/(8*r_h^2);
  result.d_h := d_h;
  result.v_mean := u;
  result.A := pi*(d_h^2)/4;

  annotation (Documentation(info="<html>
<p>Pressure loss after after&nbsp;Darcy&ndash;Weisbach, which is valid in laminar and turbulent flow regimes.</p>
<p>In order to avoid a 0^0 for Re=0 (and therefore a = 1) in the computation of lambda_aux, we add epsilon=0.01 to Re to lower bound it in a smooth way.</p>
<p>ks_input defines the pipe roughness. It can be selected from a list of materials or given directly.</p>
<p><img src=\"modelica://Thermofluidstream/Resources/Doku/ThermofluidStream.Processes.Internal.FlowResistance.laminarTurbulentPressureLoss.PNG\"/></p>
<p><br>Cheng, Nian-Sheng (2008). Formulas for friction factor in transitional regimes. In:Journal of Hydraulic Engineering134.9, pp. 1357-1362</p>
<p>Elmqvist, Hilding, Hubertus Tummescheit, and Martin Otter (2003). Object-orientedmodeling of thermo-fluid systems. In:3rd International Modelica Conference,pp. 269-286.</p>
</html>"));
end laminarTurbulentPressureLoss;
