﻿within ThermofluidStream.Processes.Internal.FlowResistance;
function laminarTurbulentPressureLossHaaland "Laminar and turbulent flow regimes pressure loss function (Haaland 1983)"
  extends Internal.FlowResistance.partialPressureLoss;

  import Modelica.Constants.pi;

  // Inputs
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

  input Real Re_laminar(unit "1") = 2000
    "Upper Reynolds number boundary for laminar flow in pipe"
    annotation(Dialog(enable=true));
  input Real Re_turbulent(unit "1") = 4000
    "Lower Reynolds number boundary for turbulent flow in pipe"
    annotation(Dialog(enable=true));
  input Real shape_factor(unit "1") = 64
    "Laminar pressure loss factor based on Hagen-Poiseuille loss"
    annotation(Dialog(enable=true));
  input Real n(unit "1") = 1 "Transition Coefficient see Documentation"
    annotation(Dialog(enable=true));

  input SI.Length ks_input(min=1e-7) = 1e-7 "Pipe roughness" annotation (Dialog(
        enable=(material == ThermofluidStream.Processes.Internal.Material.other)));

  input ThermofluidStream.Processes.Internal.Material material=
      ThermofluidStream.Processes.Internal.Material.other "Material of pipe"
    annotation (Dialog(enable=true), choices(
      choice=ThermofluidStream.Processes.Internal.Material.concrete
        "Concrete ks=5mm",
      choice=ThermofluidStream.Processes.Internal.Material.wood "Wood ks=0.5mm",
      choice=ThermofluidStream.Processes.Internal.Material.castIron
        "Cast Iron ks=0.25mm",
      choice=ThermofluidStream.Processes.Internal.Material.galvanizedIron
        "Galvanized Iron ks=0.15mm",
      choice=ThermofluidStream.Processes.Internal.Material.steel
        "Steel ks=0.059mm",
      choice=ThermofluidStream.Processes.Internal.Material.drawnPipe
        "Drawn Pipe ks=0.0015mm"));

protected
  SI.Length ks "pipe roughness";

  SI.Area area=pi*(d_h^2)/4 "Area of Pipe";
  Real relative_roughness=ks/d_h "Relative Roughness of Pipe";

  Real Re_abs "absolute value of Reynolds number";
  Real Re_abs_limited "limited absolute value of Reynolds number";

  Real friction_factor;
  SI.Pressure pressureLossLaminar "Laminar Pressure Loss";
  SI.Pressure pressureLossTurbulent "Turbulent Pressure Loss";

  constant Real eps=1e-5 "lower bound of tubulent Re to avoid devision by zero";

algorithm
  if material == ThermofluidStream.Processes.Internal.Material.concrete then
    ks := 5e-3;
  elseif material == ThermofluidStream.Processes.Internal.Material.wood then
    ks := 0.5e-3;
  elseif material == ThermofluidStream.Processes.Internal.Material.castIron then
    ks := 0.25e-3;
  elseif material == ThermofluidStream.Processes.Internal.Material.galvanizedIron then
    ks := 0.15e-3;
  elseif material == ThermofluidStream.Processes.Internal.Material.steel then
    ks := 0.059e-3;
  elseif material == ThermofluidStream.Processes.Internal.Material.drawnPipe then
    ks := 0.0015e-3;
  else
    ks := ks_input;
  end if;

  if shape == ThermofluidStream.Processes.Internal.ShapeOfResistance.circular then
    d_h := 2*r;
  elseif shape == ThermofluidStream.Processes.Internal.ShapeOfResistance.rectangle then
    d_h := 2*a*b/(a+b);
  elseif shape == ThermofluidStream.Processes.Internal.ShapeOfResistance.other then
    d_h := d_h_input;
  end if;

  // absolute Reynolds number
  Re_abs := abs(m_flow)*d_h/(area*mu);
  Re_abs_limited := max(eps, min(1, Re_abs));

  friction_factor :=
    (-1.8/n*log10((6.9/Re_abs_limited)^n + (relative_roughness/3.75)^(1.11*n)))^(-2);

  pressureLossLaminar := m_flow*mu*shape_factor*l/(2*rho*d_h^2*area);
  pressureLossTurbulent := m_flow*abs(m_flow)*friction_factor*l/(2*rho*d_h*
    area^2);

  result.dp := Utilities.Functions.blendFunction(
    y1=pressureLossLaminar,
    y2=pressureLossTurbulent,
    x1=Re_laminar,
    x2=Re_turbulent,
    x=Re_abs);

  result.A := area;
  result.d_h := d_h;
  result.v_mean := m_flow/(rho*area);

  annotation (Documentation(info="<html>
<p>Pressure loss after&nbsp;Darcy&ndash;Weisbach, which is valid in laminar and turbulent flow regimes.</p>
<p>The friction factor is based on Haaland 1983.</p>
<p>Haaland, S. E. Simple and Explicit Formulas for the Friction Factor in Turbulent Pipe Flow. Journal of Fluids Engineering 105, 89–90; 10.1115/1.3240948 (1983).</p>
</html>"));
end laminarTurbulentPressureLossHaaland;
