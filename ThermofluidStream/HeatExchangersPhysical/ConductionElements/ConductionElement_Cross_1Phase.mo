within ThermofluidStream.HeatExchangersPhysical.ConductionElements;
model ConductionElement_Cross_1Phase "Cross Flow ConductionElement for single-phase fluids"
  extends CrossIcon;
  extends Internal.PartialConductionCell;

  import Modelica.Units.SI;

  parameter SI.Length d_out = 0.012   "Outer Pipe Diameter";
  parameter SI.Length L_pipe = 0.28   "Length of flow";
  parameter SI.Length H_ch = 0.024    "Height of flow channel";
//  parameter SI.Area   A_HX = 0.010    "heat exchange area";
  parameter Integer   N_fin = 25      "Number of fins";

  final parameter SI.Area   A_HX_pipe_out = Modelica.Constants.pi * d_out * L_pipe  "heat exchange area outer pipe surface";
  final parameter SI.Area   A_HX_parallel = 2 * N_fin * H_ch * L_pipe               "heat exchange area parallel to flow";
  final parameter SI.Length L_around =   d_out * Modelica.Constants.pi / 2;
  final parameter SI.Area   A_flow = H_ch * L_pipe;
  final parameter SI.Length d_hyd = 2*A_flow / (H_ch + L_pipe)  "hydraulic diameter of rectangular flow channel";
//  final parameter SI.Volume V = A_flow * L_pipe;

  constant Real Re_exp(unit="1") = 0.8 "Reynolds-exponent for heat transfer calculation";

  constant Real cd_cyl(unit="1") = 1.2 "drag coefficient of flow around a cylinder";

  SI.Temperature T_in    "inflow temperature";
  SI.Temperature T_wall  "wall temperature";

  SI.Velocity v_mean     "mean fluid velocity";
  SI.Velocity v_around   "mean fluid velocity around pipe";

  SI.DynamicViscosity eta_in  "dynamic viscosity";

  SI.AbsolutePressure p_dyn(displayUnit="Pa")  "dynamic pressure";

  SI.ReynoldsNumber Re_pipe "Reynolds number";
  SI.ReynoldsNumber Re_rect "Reynolds number in rectangular flow channel";

  SI.PrandtlNumber Pr_F  "Prandtl Number at inlet temperature";
//  SI.PrandtlNumber Pr_W  "Prandtl Number at wall temperature";

//  SI.PecletNumber Pe "Peclet number";

  SI.NusseltNumber Nu_pipe    "Nusselt Number for heat flow to pipe";
  SI.NusseltNumber Nu_plate   "Nusselt Number for heat flow to plates";

  Real Kfac = 1  "K factor for liquid fluids";

  SI.CoefficientOfHeatTransfer alpha_pipe   "Coefficient of Heat transfer to pipe";
  SI.CoefficientOfHeatTransfer alpha_plate  "Coefficient of Heat transfer to plate";

  SI.ThermalConductance kA_pipe       "Thermal conductance fluid -> pipe wall";
  SI.ThermalConductance kA_plate      "Thermal conductance fluid -> plates";

  SI.HeatFlowRate Q_flow_pipe     "Heat Flow to pipe surface";
  SI.HeatFlowRate Q_flow_plate    "Heat Flow to plates";

  SI.EnthalpyFlowRate H_flow_in   "Enthalpy Flow at inlet";
  SI.EnthalpyFlowRate H_flow_out  "Enthalpy Flow at outlet";

  SI.SpecificHeatCapacityAtConstantPressure cp_F   "specific heat capacity of medium";

  SI.ThermalConductivity lambda_F  "Thermal Conductivity of Fluid";

  Real eps   "Hollow Space Ratio";

equation

// Inflow and Wall resp. Tin Temperature
   T_in   = Medium.temperature(inlet.state);
   T_wall = heatPort.T;

// Enthalpy Inflow
   H_flow_in = h_in * m_flow;

// Dynamic Viscosity at inflow
   eta_in = Medium.dynamicViscosity(inlet.state);

// Hollow Space Ratio
   eps = 1 - (Modelica.Constants.pi/4) * (d_out/H_ch);

// Mean Fluid Velocity in flow channel
   v_mean = inlet.m_flow/(rho*A_flow);

// Mean Fluid Velocity around pipe
   v_around = v_mean / eps;

// Reynolds Numbers
   Re_pipe = max(1.0,abs(v_around) * L_around * rho/eta_in);
   Re_rect = max(1.0,abs(v_mean)   * d_hyd    * rho/eta_in);

// Prandtl Number at inflow
   Pr_F = Medium.prandtlNumber(inlet.state);
//  Pr_W = Medium.prandtlNumber(Medium.setState_pTX(p_in, T_wall, Xi_in));

//  Pe = Re * Pr_F;

// Nusselt Number for pipes - TODO: need to be corrected
   if Re_pipe < 2300 then
     Nu_pipe  = 0.664 * Re_pipe^(1/2) * Pr_F^(1/3) * Kfac;
   else
     Nu_pipe  = 0.037 * Re_pipe^0.8 * Pr_F / (1+2.443*Re_pipe^(-0.1) * (Pr_F^(2/3)-1)) * Kfac;
   end if;

// Nusselt Number for plates
   if Re_rect < 2300 then
     Nu_plate = 0.664 * Re_rect^(1/2) * Pr_F^(1/3) * Kfac;
   else
     Nu_plate = 0.037 * Re_rect^0.8 * Pr_F / (1+2.443*Re_rect^(-0.1) * (Pr_F^(2/3)-1)) * Kfac;
   end if;

// Thermal Conductivity at inflow
   lambda_F = Medium.thermalConductivity(inlet.state);

// Heat Capacity at inflow
   cp_F = Medium.specificHeatCapacityCp(inlet.state);

  alpha_pipe  = lambda_F/d_out * Nu_pipe;
  alpha_plate = lambda_F/d_hyd * Nu_plate;

  kA_pipe  = alpha_pipe  * A_HX_pipe_out;
  kA_plate = alpha_plate * A_HX_parallel;

  Q_flow_pipe  = kA_pipe * (T_wall - T_in);
  Q_flow_plate = if noEvent(m_flow > 0) then cp_F * m_flow * (T_wall - T_in) * (1 - exp( -kA_plate / (cp_F * m_flow))) else kA_plate*(T_wall - T_in);

  Q_flow = Q_flow_pipe + Q_flow_plate;

// No Change of Mass Flow
// ----------------------
   dm_flow = 0;

// No Change of Composition
   Xi_out = Xi_in;

// Change of specific Enthalpy (see PartialConductionCell)
// -------------------------------------------------------
   h_out = h;

   H_flow_out = H_flow_in + Q_flow;

// -------------------------------------------------------------------------------
// Pressure Drop against a crosswise pipe
//--------------------------------------------------------------------------------

// Dynamic Fluid Pressure around pipe
   p_dyn = rho/2 * v_around^2;

// Pressure Drop from drag force of a cylinder
   dp = - cd_cyl * d_out/H_ch * p_dyn;

  annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
        Text(
          extent={{-68,-72},{70,-100}},
          textColor={0,0,0},
          textString="1-Phase Medium")}),
    Diagram(coordinateSystem(preserveAspectRatio=false)),
    Documentation(info="<html>
    For a single conduction element this component models
    <ul>
    <li>the temperature drop or rise of the fluid along the cross channel</li>
    <li>the heat flow between fluid and pipe wall resp. fins</li>
    <li>the pressure drop along the cross channel.</li>
    </ul>

    <p>
    <strong>Media</strong>
    <p>
    The fluid medium must keep its phase along the cross flow channel, i.e. it must be either liquid or gas.<br>
    Possible fluids are liquid water, glycol, thermal oil, non-condensing gas.

    <p>
    <strong>Modeling Approach</strong>
    <p>
    The heat exchange rate is calculated according to the Nusselt Number.<br>
    The Nusselt number for flows around pipes is calculated according to Gnielinksi<br>
    as documented with formula (3.211) in the text book from H.D. Baehr and K. Stephan: \"W&auml;rme- und Stoff&uuml;bertragung\", 7th edition.<br>
    <p>
    The total heat exchange rate is the sum of heat exchange to the pipe and to the fins.<br>
    The heat exchange rate to the fins is calulated at the beginning of the fins. Along the length of the fins an exponential convergence to the wall temperature is assumed.
</html>"));
end ConductionElement_Cross_1Phase;
