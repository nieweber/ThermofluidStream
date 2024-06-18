within ThermofluidStream.HeatExchangersPhysical.ConductionElements;
model ConductionElement_Pipe_CondensingAir "Pipe Flow ConductionElement for Air with condensing Humidity"
  extends PipeIcon;
  extends Internal.PartialConductionCell(
                                redeclare package Medium =
        ThermofluidStream.Media.myMedia.Air.MoistAir);
//  extends PartialConductionCell(redeclare package Medium =
//        ThermofluidStream.Media.myMedia.Interfaces.PartialCondensingGases);  // Moist Air

  replaceable package Medium_Water = ThermofluidStream.Media.myMedia.Water.StandardWater
    "Medium model"
    annotation (choicesAllMatching=true, Dialog(group = "Medium definitions"));

  import Modelica.Units.SI;
  import ThermofluidStream.Media.myMedia.Water.StandardWater;
  import ThermofluidStream.Media.myMedia.Air.DryAirNasa;

  parameter SI.Length d_hyd = 1   "Hydraulic Diameter";
  parameter SI.Length L_flow = 1  "Length of flow";

  parameter SI.Duration dt = 1  "reaction time";

  final parameter SI.Area A_flow = Modelica.Constants.pi * (d_hyd/2)^2;
  final parameter SI.Area A_inside =  Modelica.Constants.pi * d_hyd * L_flow;

//  constant SI.Duration dt_cond = 0.1  "time constant of condensation";

  constant Real Re_exp(unit="1") = 0.8 "Reynolds-exponent for heat transfer calculation";

  constant SI.SpecificHeatCapacity R_da =  Medium.dryair.R_s  "Ideal gas constant of dry air";
  constant SI.SpecificHeatCapacity R_vap = Medium.steam.R_s   "Ideal gas constant of water vapor";

  constant Real Le_vap = 0.87  "Lewis number for water vapor";

  constant SI.ReynoldsNumber Re_crit_low =  1500  "lower limit of transition phase laminar to turbulent";
  constant SI.ReynoldsNumber Re_crit_high = 2500  "upper limit of transition phase laminar to turbulent";

  SI.Density rho_in   "2-phase density at inlet";
  SI.Density rho_gas  "density in liquid phase at inlet pressure";
  SI.Density rho_liq  "density in liquid phase at inlet pressure";

  SI.AbsolutePressure p_vap_in        "partial vapor pressure at inlet";
  SI.AbsolutePressure p_vap_sat_wall  "vapor saturation pressure at wall temperature";
  SI.AbsolutePressure p_vap_sat_0     "vapor saturation pressure at condensation temperature iteration start value";
  SI.AbsolutePressure p_vap_sat_1     "vapor saturation pressure at condensation temperature in 1st iteration step";
  SI.AbsolutePressure p_vap_sat_2     "vapor saturation pressure at condensation temperature in 2nd iteration step";

  SI.Temperature T_in        "inflow temperature";
  SI.Temperature T_wall      "wall temperature";
  SI.Temperature T_sat       "saturation temperature of condensing gas";
  SI.Temperature T_cond_0    "condensation temperature start value";
  SI.Temperature T_cond_1    "condensation temperature 1st iteration";
  SI.Temperature T_cond_2    "condensation temperature 2nd iteration";
  SI.Temperature T_cond_mean "condensation temperature mean value";
  SI.Temperature T_water     "Mean Temperature of condensed water";

  Real phi_in  "relative humdity inflow";

  SI.SpecificEnthalpy dh_vap  "enthalpy of vaporization";

  SI.Velocity v_mean     "mean fluid velocity";

  SI.DynamicViscosity mu_gas_in   "dynamic viscosity of gas at inlet";
  SI.DynamicViscosity mu_liq      "dynamic viscosity of liquid water at wall temperature";

  SI.ReynoldsNumber Re "Reynolds number";

  SI.PrandtlNumber Pr_gas  "Prandtl Number at inlet temperature";

  SI.PecletNumber Pe "Peclet number";

  SI.NusseltNumber Nu            "Nusselt Number";

  Real zeta   "Pipe friction zeta value";

  SI.CoefficientOfHeatTransfer beta_gas  "Mass Transfer Coefficient";

  Real Phi_Ack   "Coefficient of Ackermann Correction";
  Real zeta_Ack  "Ackermann Correction";

  SI.CoefficientOfHeatTransfer alpha_mean          "Mean Coefficient of Heat transfer";
  SI.CoefficientOfHeatTransfer alpha_gas           "Heat Transfer Coefficient of gas to wall and film";
  SI.CoefficientOfHeatTransfer alpha_gas_Ack       "Heat Transfer Coefficient of gas to wall and film with Ackermann Correction";
  SI.CoefficientOfHeatTransfer alpha_cond          "Heat Transfer Coefficient by condensation";
  SI.CoefficientOfHeatTransfer alpha_cond_1        "Heat Transfer Coefficient by condensation in 1st iteration step";
  SI.CoefficientOfHeatTransfer alpha_cond_2        "Heat Transfer Coefficient by condensation in 2nd iteration step";

  SI.ThermalConductivity lambda_gas  "Thermal Conductivity of inflow gas";
  SI.ThermalConductivity lambda_liq  "Thermal Conductivity in liquid state";

  SI.SpecificHeatCapacityAtConstantPressure cp_gas_in  "specific heat capacity of fluid at inlet";
  SI.SpecificHeatCapacityAtConstantPressure cp_air_in  "specific heat capacity of dry air at inlet";
  SI.SpecificHeatCapacityAtConstantPressure cp_vap_in  "specific heat capacity of water vapor at inlet";

  Real X_humid_in      "inflow water partition: vapor + liquid water";
  Real X_vap_sat_in    "maximum vapor partition at inflow";
  Real X_vap_in        "inflow vapor partition";
  Real X_liq_in        "inflow liquid water partition";
  Real X_vap_sat_wall  "minimum vapor partition at wall temperature";
  Real X_humid_out     "outflow water partition: vapor + liquid water";
  Real X_vap_sat_out   "maximum vapor partition at outflow";
  Real X_vap_out       "vapor partition at outflow";
  Real X_liq_out       "outflow liquid water partition";
  Real X_vap_cond      "vapor partition at condensation area";

  SI.Length delta_film_mean_1     "Mean Thickness of condensed film in pipe 1st iteration step";
  SI.Length delta_film_mean_2     "Mean Thickness of condensed film in pipe 2nd iteration step";

  SI.MassFlowRate m_flow_cond_0            "condensing mass flow iteration start value";
  SI.MassFlowRate m_flow_cond_1            "condensing mass flow 1st iteration step";
  SI.MassFlowRate m_flow_cond_2            "condensing mass flow 2nd iteration step";
  SI.MassFlowRate m_flow_cond_cell         "average condensing mass flow in cell";
//  SI.MassFlowRate m_flow_cond_equlib       "mass flow rate of flowing condensate in equilibrium (Nusselt Formula)";

  SI.MassFlowRate m_flow_out          "outlet mass flow";
  SI.MassFlowRate m_flow_humid_out    "humid outlet mass flow (vapor + liquid";

  SI.HeatFlowRate Q_flow2           "heat flow from gas inflow to wall";
  SI.HeatFlowRate Q_flow_cond       "heat flow of condensation";
  SI.HeatFlowRate Q_flow_gas_film   "heat flow from gas to condensate film";
  SI.HeatFlowRate Q_flow_gas_film2  "heat flow from gas inflow to condensate film";
  SI.HeatFlowRate Q_flow_gas_wall   "heat flow from gas to wall";
  SI.HeatFlowRate Q_flow_gas_wall2  "heat flow from gas inflow to wall";

  SI.EnthalpyFlowRate H_flow_in        "Enthalpy Flow at inlet";
  SI.EnthalpyFlowRate H_flow_out       "Enthalpy Flow at outlet";
  SI.EnthalpyFlowRate H_flow_water     "Enthalpy Flow of condensed water";

  SI.SpecificEnthalpy h_water          "Specific Enthalpy of condensed water";

  Real m_flow_density(unit="kg/(s.m2)")      "mass flow density";

  Real X_lam   "laminar portion";
  Real X_turb  "turbulent portion";

  Real pipeFric_gas    "Pipe friction lambda value in gas phase";

  Real pipeFric_lam    "Pipe friction lambda value in laminar flow";
  Real pipeFric_turb   "Pipe friction lambda value in turbulent flow";

  Real dpdz_fric(unit="Pa/m")   "pressure gradient by friction";

  ThermofluidStream.Topology.JunctionT1 junctionT1_1(redeclare package Medium = Medium_Water)                                        annotation (Placement(transformation(extent={{10,-70},{-10,-50}})));
  ThermofluidStream.Boundaries.Source source(
    redeclare package Medium = Medium_Water,
    setEnthalpy=false,
    pressureFromInput=true,
    temperatureFromInput=true,
    enthalpyFromInput=false)                                                                                                 annotation (Placement(transformation(
        extent={{-10,10},{10,-10}},
        rotation=-90,
        origin={0,20})));
  Modelica.Blocks.Sources.RealExpression MassFlow_condensed(y=-dm_flow)    annotation (Placement(transformation(extent={{-90,30},{-70,50}})));
  Modelica.Blocks.Sources.RealExpression Temperature_Water_condensed(y=T_water) annotation (Placement(transformation(extent={{34,30},{14,50}})));
  ThermofluidStream.Sensors.SingleFlowSensor condensedWaterFlowSensor(
    redeclare package Medium = Medium_Water,
    quantity=ThermofluidStream.Sensors.Internal.Types.MassFlowQuantities.m_flow_kgps,
    outputValue=true) annotation (Placement(transformation(
        extent={{-10,10},{10,-10}},
        rotation=-90,
        origin={-6,-20})));
  HXutilities.Controller.PID_varLimits condensedWater_PressureController(k=1000, xi_start=1.0E+05) annotation (Placement(transformation(extent={{-60,30},{-40,50}})));
  ThermofluidStream.Interfaces.Inlet inlet_water(redeclare package Medium = Medium_Water) annotation (Placement(transformation(extent={{-120,-80},{-80,-40}}), iconTransformation(extent={{-120,-80},{-80,-40}})));
  ThermofluidStream.Interfaces.Outlet outlet_water(redeclare package Medium = Medium_Water) annotation (Placement(transformation(extent={{80,-80},{120,-40}}), iconTransformation(extent={{80,-80},{120,-40}})));

initial equation
   dm_flow = 0.0;

equation

  T_in   = Medium.temperature(inlet.state);
  T_wall = heatPort.T;

// Enthalpy Inflow
   H_flow_in = h_in * m_flow;

  phi_in = Medium.relativeHumidity(inlet.state);

  X_humid_in   = inlet.state.X[Medium.Water];
  X_vap_sat_in = Medium.Xsaturation(inlet.state);
  X_vap_in     = min(X_humid_in, X_vap_sat_in);
  X_liq_in     = X_humid_in - X_vap_in;

  X_vap_sat_wall = Medium.Xsaturation(Medium.setState_pTX(p_in, T_wall));

  // The function specificHeatCapacityCp does not yield correct results if X_liq_in > 0
  //cp_gas_in = Medium.specificHeatCapacityCp(inlet.state);

  cp_air_in = DryAirNasa.specificHeatCapacityCp(inlet.state);
  cp_vap_in = ThermofluidStream.Media.myMedia.IdealGases.SingleGases.H2O.specificHeatCapacityCp(inlet.state);

  cp_gas_in = (cp_air_in * (1 - X_vap_in) + cp_vap_in * X_vap_in) / (1-X_liq_in);

  rho_in = Medium.density(inlet.state);

  lambda_gas = Medium.thermalConductivity(inlet.state);

  v_mean = inlet.m_flow/(rho_in*A_flow);

  mu_gas_in = Medium.dynamicViscosity(inlet.state);

  Re = max(1.0,abs(v_mean) * d_hyd * rho_in/mu_gas_in);

  Pr_gas = mu_gas_in * cp_gas_in / lambda_gas;

// Pr_gas = Medium.prandtlNumber(inlet.state);
// not usable since cp is not calculated correctly

//------------------------------------------------------------------------
// Heat Exchange without Condensation
//------------------------------------------------------------------------

  if noEvent(Re < 2300) then   // Laminar Flow

    zeta = 0;
    Pe   = Re * Pr_gas;

    Nu = (49.37 + (1.615*(Pe*d_hyd/L_flow)^(1/3) -0.7)^3)^(1/3);

  else   // Turbulent Flow

    Pe = 0;

    zeta = (0.79 * log(Re) - 1.64)^(-2);
    Nu = zeta/8*(Re-1000)*Pr_gas / (1+12.7*(Pr_gas^(2/3)-1)*sqrt(zeta/8)) * (1+(d_hyd/L_flow)^(2/3));

  end if;

  alpha_gas = lambda_gas/d_hyd * Nu;

  Q_flow_gas_wall = if noEvent(m_flow > 0) then cp_gas_in * m_flow * (T_wall - T_in) * (1 - exp( -(alpha_gas * A_inside) / (cp_gas_in * m_flow))) else 0;

  Q_flow_gas_wall2 = alpha_gas * A_inside * (T_wall - T_in);

//------------------------------------------------------------------------
// Check if condensation can take place
//------------------------------------------------------------------------

  p_vap_in = X_vap_in * R_vap / (R_da * (1-X_vap_in) + X_vap_in * R_vap) * p_in;

  p_vap_sat_wall = Medium.saturationPressureLiquid(T_wall);

  if noEvent(p_vap_in > p_vap_sat_wall + 1/1000) then

//------------------------------------------------------------------------
// Material Constants under current conditions
//------------------------------------------------------------------------

    // Gas Density
    rho_gas = ThermofluidStream.Media.myMedia.Air.DryAirNasa.density_pT(p_in, T_wall);

    // Material Coefficients of liquid water at wall conditions
    rho_liq    = StandardWater.density_pT(p_in, T_wall);
    lambda_liq = StandardWater.thermalConductivity(StandardWater.setState_pTX(p_in, T_wall));
    mu_liq     = StandardWater.dynamicViscosity(StandardWater.setState_pTX(p_in, T_wall));

    // Enthalpy of Condensation
    dh_vap = Medium.enthalpyOfVaporization(T_wall);

    beta_gas = alpha_gas / (cp_gas_in * rho_in) * Le_vap^(-2/3);

    // Saturation Temperature at inflow Pressure
    T_sat = Medium.saturationTemperature(p_vap_in);

//------------------------------------------------------------------------
// Vapor Transfer to Condensation Area
//------------------------------------------------------------------------

    T_cond_0 = T_wall;

    p_vap_sat_0 = p_vap_sat_wall;

    m_flow_cond_0 = alpha_gas / cp_gas_in * Le_vap^(-2/3) * A_inside * (p_vap_in - p_vap_sat_0) / (p_in - p_vap_sat_0);

    // 1st Iteration -------------------------------------------------------

    // Ackermann Correction
    Phi_Ack       = abs(m_flow_cond_0) * cp_gas_in / (A_inside * alpha_gas);
    zeta_Ack      = if noEvent(Phi_Ack > 0) then - Phi_Ack / (exp(-Phi_Ack) - 1) else 1;
    alpha_gas_Ack = zeta_Ack * alpha_gas;

    // Theoretic Mean Film Thickness
    delta_film_mean_1 = ((3 * mu_liq * m_flow_cond_0) / (Modelica.Constants.g_n * rho_liq * (rho_liq-rho_gas) * L_flow))^(1/3);

    alpha_cond_1 = lambda_liq / delta_film_mean_1;

    T_cond_1 = min(T_sat, max(T_wall, (dh_vap * m_flow_cond_0 / A_inside + alpha_gas_Ack * T_in + alpha_cond_1 * T_wall) / (alpha_gas_Ack + alpha_cond_1)));

    p_vap_sat_1 = min(p_vap_in, Medium.saturationPressureLiquid(T_cond_1));

    m_flow_cond_1 = alpha_gas_Ack / cp_gas_in * Le_vap^(-2/3) * A_inside * (p_vap_in - p_vap_sat_1) / (p_in - p_vap_sat_1);

    // 2nd Iteration -------------------------------------------------------
    if noEvent(m_flow_cond_1 > 0 and T_cond_1 - T_cond_0 > 0.1) then

      // Theoretic Mean Film Thickness on a vertical plate
      delta_film_mean_2 = ((3 * mu_liq * m_flow_cond_1) / (Modelica.Constants.g_n * rho_liq * (rho_liq-rho_gas) * L_flow))^(1/3);

      alpha_cond_2 = lambda_liq / delta_film_mean_2;

      T_cond_2 = min(T_sat, max(T_wall, (dh_vap * m_flow_cond_1 / A_inside + alpha_gas_Ack * T_in + alpha_cond_2 * T_wall) / (alpha_gas_Ack + alpha_cond_2)));

      p_vap_sat_2 = min(p_vap_in, Medium.saturationPressureLiquid(T_cond_2));

      m_flow_cond_2 = alpha_gas_Ack / cp_gas_in * Le_vap^(-2/3) * A_inside * (p_vap_in - p_vap_sat_2) / (p_in - p_vap_sat_2);

    else
      delta_film_mean_2 = delta_film_mean_1;
      alpha_cond_2      = alpha_cond_1;
      T_cond_2          = T_cond_1;
      p_vap_sat_2       = p_vap_sat_1;
      m_flow_cond_2     = m_flow_cond_1;
    end if;

    X_vap_cond = R_da/R_vap * p_vap_sat_2/p_in;

    m_flow_cond_cell = if noEvent(m_flow > 0) then m_flow * (X_vap_in - X_vap_cond) * (1 - exp(-(alpha_gas_Ack * R_vap * A_inside) / (cp_gas_in * Le_vap^(2/3) * R_da * m_flow))) else 0;

    T_cond_mean = T_cond_2;

    Q_flow_cond = - dh_vap *  m_flow_cond_cell;

    Q_flow_gas_film = if noEvent(m_flow > 0) then cp_gas_in * m_flow * (T_cond_mean  - T_in) * (1 - exp( -(alpha_gas_Ack * A_inside) / (cp_gas_in * m_flow))) else 0;

    Q_flow_gas_film2 = alpha_gas_Ack * (T_cond_mean - T_in) * A_inside;

    Q_flow = Q_flow_cond + Q_flow_gas_film;

    Q_flow2 = Q_flow_cond + Q_flow_gas_film2;

    alpha_cond = -Q_flow_cond / (T_cond_mean - T_wall) / A_inside;

    // Estimated temperature of condensed water
    T_water = if noEvent(m_flow_cond_cell > 0) then (T_cond_mean + T_wall)/2 else T_wall;

    // Specific enthalpy of condensed water
    h_water = StandardWater.specificEnthalpy_pT(p_in, T_water);

    // Enthalpy of outflow gas
    //H_flow_out = H_flow_in + Q_flow_gas_film;

    // Enthalpy Flow of water
    H_flow_water = h_water * m_flow_cond_cell;

    // Maximum condensing mass flow (equilibrium conditions)
//    m_flow_cond_equlib = 2.28714 * L_flow * ((d_hyd * (T_sat - T_wall) * lambda_liq) / dh_vap)^(3/4) * ((Modelica.Constants.g_n * rho_liq * (rho_liq-rho_gas)) / mu_liq)^(1/4);

//------------------------------------------------------------------------
  else    //  Heat Flow of Humid Air without Condensation
//------------------------------------------------------------------------

    T_cond_0    = 0;
    T_cond_1    = 0;
    T_cond_2    = 0;
    T_cond_mean = 0;

    p_vap_sat_0 = 0;
    p_vap_sat_1 = 0;
    p_vap_sat_2 = 0;

    m_flow_cond_0 = 0;
    m_flow_cond_1 = 0;
    m_flow_cond_2 = 0;
    m_flow_cond_cell = 0;

    delta_film_mean_1 = 0;
    delta_film_mean_2 = 0;

    alpha_cond_1 = 0;
    alpha_cond_2 = 0;

    X_vap_cond = 0;

    beta_gas = 0;

    Phi_Ack  = 0;
    zeta_Ack = 0;
    alpha_gas_Ack = alpha_gas;

    rho_gas = 0;
    rho_liq = 0;
    lambda_liq = 0;
    mu_liq = 0;
    dh_vap = 0;
    T_sat = 0;

    Q_flow  = Q_flow_gas_wall;
    Q_flow2 = Q_flow_gas_wall2;
    Q_flow_cond = 0;
    Q_flow_gas_film = 0;
    Q_flow_gas_film2 = 0;

//    m_flow_cond_equlib = 0;

    T_water = 290.0;

    alpha_cond = 0;

    h_water = 0.0;

    //H_flow_out   = H_flow_in + Q_flow_gas_wall;
    H_flow_water = 0.0;

  end if;

//------------------------------------------------------------------------
// Final Air and Water Mass Flow
//------------------------------------------------------------------------

   der(dm_flow) = (-m_flow_cond_cell - dm_flow) / dt;

   m_flow_out = m_flow + dm_flow;

   m_flow_humid_out = X_humid_in * m_flow + dm_flow;   // TODO m_flow = m_flow_in

//------------------------------------------------------------------------
// Final Outflow Humdity
//------------------------------------------------------------------------

   X_humid_out = if noEvent(m_flow_out > 0) then m_flow_humid_out / m_flow_out else X_humid_in;

   Xi_out = {X_humid_out};

// Outflow water and vapor Partition
// ---------------------------------
   X_vap_sat_out = Medium.Xsaturation(outlet.state);
   X_vap_out     = min(X_humid_out, X_vap_sat_out);
   X_liq_out     = X_humid_out - X_vap_out;

//------------------------------------------------------------------------
// Final Outflow Heat
//------------------------------------------------------------------------

   H_flow_out = H_flow_in + Q_flow - H_flow_water;

// Change of specific Enthalpy (see PartialConductionCell)
// -------------------------------------------------------
   h_out = if noEvent(m_flow_out > 0 and H_flow_out >0) then H_flow_out / m_flow_out else h_in;

   //Q_flow = -(H_flow_in - H_flow_out - H_flow_water);

// Mean Heat Transfer Coefficient
// ------------------------------
   alpha_mean = -Q_flow / (T_in - T_wall) / A_inside;

// -------------------------------------------------------------------------------
// Pressure Drop
//--------------------------------------------------------------------------------

// Mass Flow Density
// -----------------
   m_flow_density = m_flow/A_flow;

// Laminar/turbulent partitions
// ----------------------------
   (X_lam, X_turb) =Internal.portion(Re_crit_low, Re_crit_high, Re);

// Friction Factors in laminar and turbulent flow
// ----------------------------------------------
   pipeFric_lam  = 64/Re;
   pipeFric_turb = 0.3164/(Re^0.25);   //smooth pipe

// Combined Friction Factor
   pipeFric_gas = X_lam * pipeFric_lam + X_turb * pipeFric_turb;

// 1 Phase Flow Friction Pressure Drop Gradient
// --------------------------------------------
   dpdz_fric = -pipeFric_gas/d_hyd / (2*rho) * m_flow_density * abs(m_flow_density);

// Pressure drop proportional to length of pipe
// --------------------------------------------
   dp = dpdz_fric * L_flow;

  connect(junctionT1_1.outlet,outlet_water)  annotation (Line(
      points={{10,-60},{100,-60}},
      color={28,108,200},
      thickness=0.5));
  connect(inlet_water,junctionT1_1. inletB) annotation (Line(
      points={{-100,-60},{-40,-60},{-40,-88},{0,-88},{0,-70}},
      color={28,108,200},
      thickness=0.5));
  connect(Temperature_Water_condensed.y,source. T0_var) annotation (Line(points={{13,40},{0,40},{0,22}}, color={0,0,127}));
  connect(source.outlet,condensedWaterFlowSensor. inlet) annotation (Line(
      points={{0,10},{0,-10}},
      color={28,108,200},
      thickness=0.5));
  connect(condensedWaterFlowSensor.outlet,junctionT1_1. inletA) annotation (Line(
      points={{0,-30},{0,-50}},
      color={28,108,200},
      thickness=0.5));
  connect(MassFlow_condensed.y,condensedWater_PressureController. u_s) annotation (Line(points={{-69,40},{-61,40}}, color={0,0,127}));
  connect(condensedWater_PressureController.u_m,condensedWaterFlowSensor. value_out) annotation (Line(points={{-50,29},{-50,-40},{-6,-40},{-6,-30}}, color={0,0,127}));
  connect(condensedWater_PressureController.y,source. p0_var) annotation (Line(points={{-39,40},{-6,40},{-6,22}}, color={0,0,127}));
  annotation (
    Icon(coordinateSystem(preserveAspectRatio=false), graphics={Text(
          extent={{-70,-66},{70,-100}},
          textColor={0,0,0},
          textString="Condensing Air")}),
    Diagram(coordinateSystem(preserveAspectRatio=false)),
    Documentation(info="<html>
    For a single conduction element this component models
    <ul>
    <li>the temperature drop or rise of the fluid along the pipe channel</li>
    <li>the heat flow between fluid and pipe wall</li>
    <li>the pressure drop along the pipe channel.</li>
    </ul>
    Condensation of Water Vapor to Liquid Water is supported.

    <p>
    <strong>Media</strong>
    <p>
    The fluid medium is moist air. The humidity in the air may condense. The condensed water is treated as a separate water flow.

    <p>
    <strong>Modeling Approach</strong>
    <p>
    <em>Non-condensing Conditions</em>
    <p>
    The heat exchange rate is calculated according to the Nusselt Number.<br>
    The Nusselt number for flows in pipes is calculated
    <ul>
    <li>in the case of laminar flow according to schweizer-fn.de</li>
    <li>in the case of turbulent flow according to Gnielinksi as documented with formula (3.260)</li>
    as documented with formula (3.211) in the text book from H.D. Baehr and K. Stephan: \"W&auml;rme- und Stoff&uuml;bertragung\", 7th edition.</li>
    </ul>
    <p>
    The heat exchange rate to the pipe wall is calulated at the beginning of the pipe. Along the length of the pipe an exponential convergence to the wall temperature is assumed.
    <p>
    <em>Condensation</em>
    <p>
    Condensation can take place when the vapor pressure of the air exceeds the saturation vapor pressure at wall temperature.<br>
    The condensation mass flow is mainly determined by the transport of the distributed vapor to the pipe surface.
    This is a diffusion process which is driven by the vapor pressure difference as expressed in formula (4.31)<br>
    The condensed water forms a water film on the pipe surface, which affects the heat exchange rate as described in formula (4.30).<br>
    All the effects are described in chapter 4.14 \"Einfluss nicht kondensierbarer Gase\" in the named text book.
</html>"),
    __Dymola_Commands);
end ConductionElement_Pipe_CondensingAir;
