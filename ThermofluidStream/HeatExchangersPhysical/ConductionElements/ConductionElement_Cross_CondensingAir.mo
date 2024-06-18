within ThermofluidStream.HeatExchangersPhysical.ConductionElements;
model ConductionElement_Cross_CondensingAir "Cross Flow ConductionElement for Air with condensing Humidity"
  extends CrossIcon;
  extends Internal.PartialConductionCell(
                                initM_flow=ThermofluidStream.Utilities.Types.InitializationMethods.state,
                                redeclare package Medium = ThermofluidStream.Media.myMedia.Air.MoistAir); // constrainedby ThermofluidStream.Media.myMedia.Interfaces.PartialCondensingGas);
//        ThermofluidStream.Media.myMedia.Air.MoistAir);
//  extends SISOFlowCell_water(final clip_p_out_water=false);

  replaceable package Medium_Water = ThermofluidStream.Media.myMedia.Water.StandardWater
    "Medium model"
    annotation (choicesAllMatching=true, Dialog(group = "Medium definitions"));

  import Modelica.Units.SI;
  import ThermofluidStream.Media.myMedia.Water.StandardWater;
  import ThermofluidStream.Media.myMedia.Air.DryAirNasa;

  parameter SI.Length d_out = 0.012   "Outer Pipe Diameter";
  parameter SI.Length L_pipe = 0.28   "Length of flow";
  parameter SI.Length H_ch = 0.024    "Height of flow channel";
  parameter SI.Length L_ch = 0.024    "Length of flow channel including fins";
  parameter SI.Area   A_HX = 0.010    "heat exchange area";
//  parameter SI.Area   A_HX_parallel   "heat exchange area parallel to flow";
  parameter Integer   N_fin = 12      "number of vertical fins along pipe";

  parameter SI.Duration dt = 1  "reaction time";

  final parameter SI.Area   A_HX_pipe_out = Modelica.Constants.pi * d_out * L_pipe  "heat exchange area outer pipe surface";
  final parameter SI.Length L_around =   d_out * Modelica.Constants.pi / 2;
  final parameter SI.Area   A_flow = H_ch * L_pipe;
  final parameter SI.Length d_hyd = 2*A_flow / (H_ch + L_pipe)  "hydraulic diameter of rectangular flow channel";
  final parameter SI.Area   A_HX_parallel = 2 * N_fin * H_ch * L_pipe   "heat exchange area parallel to flow";
//  final parameter SI.Volume V = A_flow * L_pipe;

  constant Real Re_exp(unit="1") = 0.8 "Reynolds-exponent for heat transfer calculation";

  constant Real cd_cyl(unit="1") = 1.2 "drag coefficient of flow around a cylinder";

  SI.Temperature T_gas_in    "inflow humid air (gas) temperature";
  SI.Temperature T_wall      "wall temperature";

  SI.Density rho_in           "humid gas density at inlet";

  SI.DynamicViscosity eta_in  "dynamic viscosity";

  SI.Velocity v_mean     "mean fluid velocity";
  SI.Velocity v_around   "mean fluid velocity around pipe";

  SI.AbsolutePressure p_dyn(displayUnit="Pa")  "dynamic pressure";

  SI.ReynoldsNumber Re_pipe "Reynolds number";
  SI.ReynoldsNumber Re_rect "Reynolds number in rectangular flow channel";

  SI.PrandtlNumber Pr_F  "Prandtl Number at inlet temperature";
//  SI.PrandtlNumber Pr_W  "Prandtl Number at wall temperature";

//  SI.PecletNumber Pe "Peclet number";

  SI.NusseltNumber Nu_pipe    "Nusselt Number for heat flow to pipe";
  SI.NusseltNumber Nu_plate   "Nusselt Number for heat flow to plates";

  Real Kfac = 1  "K factor for liquid fluids";

  SI.CoefficientOfHeatTransfer alpha_air_pipe   "Coefficient of Heat transfer to pipe";
  SI.CoefficientOfHeatTransfer alpha_air_plate  "Coefficient of Heat transfer to plate";

  SI.ThermalConductance kA_pipe       "Thermal conductance fluid -> pipe wall";
  SI.ThermalConductance kA_plate      "Thermal conductance fluid -> plates";

  SI.HeatFlowRate Q_flow_pipe     "Heat Flow to pipe surface";
  SI.HeatFlowRate Q_flow_plate    "Heat Flow to plates";

  SI.SpecificHeatCapacityAtConstantPressure cp_F   "specific heat capacity of medium";

  SI.ThermalConductivity lambda_F  "Thermal Conductivity of Fluid";

  Real eps   "Hollow Space Ratio";

  // ----------------- Quantities for water condensation -----------------

  constant SI.SpecificHeatCapacity R_da =  Medium.dryair.R_s  "Ideal gas constant of dry air";
  constant SI.SpecificHeatCapacity R_vap = Medium.steam.R_s   "Ideal gas constant of water vapor";

  constant Real Le_vap = 0.87  "Lewis number for water vapor";

  SI.AbsolutePressure p_vap_in        "partial vapor pressure at inlet";
  SI.AbsolutePressure p_vap_sat_wall  "vapor saturation pressure at wall temperature";
  SI.AbsolutePressure p_vap_sat_0     "vapor saturation pressure at condensation temperature iteration start value";
  SI.AbsolutePressure p_vap_sat_pipe_1      "vapor saturation pressure on pipe at condensation temperature in 1st iteration step";
  SI.AbsolutePressure p_vap_sat_plate_1     "vapor saturation pressure on plate at condensation temperature in 1st iteration step";
  SI.AbsolutePressure p_vap_sat_pipe      "vapor saturation pressure on pipe at condensation temperature in final iteration step";
  SI.AbsolutePressure p_vap_sat_plate     "vapor saturation pressure on plate at condensation temperature in final iteration step";
//  SI.AbsolutePressure p_vap_sat_2     "vapor saturation pressure at condensation temperature in 2nd iteration step";

  SI.Temperature T_sat           "saturation temperature of condensing gas";
  SI.Temperature T_cond_0        "condensation temperature start value";
  SI.Temperature T_cond_pipe_1   "condensation temperature on plate 1st iteration";
  SI.Temperature T_cond_plate_1  "condensation temperature on plate 1st iteration";
  SI.Temperature T_cond_pipe     "condensation temperature on plate final iteration";
  SI.Temperature T_cond_plate    "condensation temperature on plate final iteration";
//  SI.Temperature T_cond_2  "condensation temperature 2nd iteration";
//  SI.Temperature T_cond_mean  "condensation temperature mean value";
  SI.Temperature T_water         "Mean Temperature of condensed water";

  SI.SpecificEnthalpy h_water_pipe    "specific enthalpy of condensed water on pipe";
  SI.SpecificEnthalpy h_water_plate   "specific enthalpy of condensed water on plate";
//  SI.SpecificEnthalpy h_out_water    "specific enthalpy of condensed water";

  SI.EnthalpyFlowRate H_flow_water_pipe    "Enthalpy Flow of condensed water on pipe";
  SI.EnthalpyFlowRate H_flow_water_plate   "Enthalpy Flow of condensed water on plate";
  SI.EnthalpyFlowRate H_flow_water         "Total Enthalpy Flow of condensed water";
  SI.EnthalpyFlowRate H_flow_in           "Enthalpy Flow at inlet";
  SI.EnthalpyFlowRate H_flow_out          "Enthalpy Flow at outlet";

  SI.MassFlowRate m_flow_vap_in   "inlet vapor flow";
//  SI.MassFlowRate m_flow_vap_out  "outlet vapor flow";

  SI.MassFlowRate m_flow_cond                "total condensing mass flow";
  SI.MassFlowRate m_flow_cond_pipe_0         "condensing mass flow to outer pipe, iteration start value";
  SI.MassFlowRate m_flow_cond_pipe_1         "condensing mass flow to outer pipe, 1st iteration step";
  SI.MassFlowRate m_flow_cond_pipe           "condensing mass flow to outer pipe, final iteration step";
//  SI.MassFlowRate m_flow_cond_2            "condensing mass flow 2nd iteration step";
//  SI.MassFlowRate m_flow_cond_cell         "average condensing mass flow in cell";

  SI.MassFlowRate m_flow_cond_plate_0       "condensing mass flow to plates, iteration start value";
  SI.MassFlowRate m_flow_cond_plate_1       "condensing mass flow to vertical fin 1st iteration step";
  SI.MassFlowRate m_flow_cond_plate         "condensing mass flow to vertical fin final iteration step";

  SI.MassFlowRate m_flow_out          "outlet mass flow";
  SI.MassFlowRate m_flow_humid_out    "humid outlet mass flow (vapor + liquid";

  SI.HeatFlowRate Q_flow_cond_pipe         "heat flow from condensed water to wall of pipe";
  SI.HeatFlowRate Q_flow_cond_plate        "heat flow from condensed water to wall of plate";
  SI.HeatFlowRate Q_flow_cond              "Total heat flow from condensed water to wall";

  Real X_humid_in(   unit="kg/kg", displayUnit="g/kg")    "inflow water partition: vapor + liquid water";
  Real X_vap_sat_in( unit="kg/kg", displayUnit="g/kg")    "maximum vapor partition at inflow";
  Real X_vap_in(     unit="kg/kg", displayUnit="g/kg")    "inflow vapor partition";
  Real X_liq_in(     unit="kg/kg", displayUnit="g/kg")    "inflow liquid water partition";
//  Real X_vap_sat_wall( unit="kg/kg", displayUnit="g/kg")  "minimum vapor partition at wall temperature";
  Real X_humid_out(    unit="kg/kg", displayUnit="g/kg")  "outflow water partition: vapor + liquid water";
//  Real X_vap_sat_out ( unit="kg/kg", displayUnit="g/kg")  "maximum vapor partition at outflow";
//  Real X_vap_out     ( unit="kg/kg", displayUnit="g/kg")  "vapor partition at outflow";
//  Real X_liq_out     ( unit="kg/kg", displayUnit="g/kg")  "outflow liquid water partition";
//  Real X_vap_cond    ( unit="kg/kg", displayUnit="g/kg")  "vapor partition at condensation area";

  SI.Density rho_gas_wall  "density in gas phase at inlet pressure and wall temperature";
  SI.Density rho_liq_wall  "density in liquid phase at inlet pressure and wall temperature";

  SI.ThermalConductivity lambda_liq  "Thermal Conductivity in liquid water at wall temperature";

  SI.DynamicViscosity eta_liq      "dynamic viscosity of liquid water at wall temperature";

  SI.SpecificEnthalpy dh_vap  "enthalpy of vaporization at wall temperature";

  SI.CoefficientOfHeatTransfer beta_air_pipe   "Mass Transfer Coefficient  TODO: Korrektur der Einheit";
  SI.CoefficientOfHeatTransfer beta_air_plate  "Mass Transfer Coefficient  TODO: Korrektur der Einheit";

  SI.SpecificHeatCapacityAtConstantPressure cp_gas_in  "specific heat capacity of fluid at inlet";
  SI.SpecificHeatCapacityAtConstantPressure cp_air_in  "specific heat capacity of dry air at inlet";
  SI.SpecificHeatCapacityAtConstantPressure cp_vap_in  "specific heat capacity of water vapor at inlet";

  Real Phi_Ack_pipe    "Phi Coefficient of Ackermann Correction (pipe)";
  Real Phi_Ack_plate   "Phi Coefficient of Ackermann Correction (plate)";

  Real zeta_Ack_pipe   "Zeta Coefficient of Ackermann Correction (pipe)";
  Real zeta_Ack_plate  "Zeta Coefficient of Ackermann Correction (plate)";

  SI.CoefficientOfHeatTransfer alpha_air_pipe_Ack     "Heat Transfer Coefficient of air to pipe wall and film with Ackermann Correction";
  SI.CoefficientOfHeatTransfer alpha_air_plate_Ack    "Heat Transfer Coefficient of air to plate wall and film with Ackermann Correction";
//  SI.CoefficientOfHeatTransfer alpha_liq_vert_0       "Heat Transfer Coefficient of liquid water on a vertical plate start value";
//  SI.CoefficientOfHeatTransfer alpha_liq_vert_1       "Heat Transfer Coefficient of liquid water on a vertical plate 1st iteration";

  SI.CoefficientOfHeatTransfer alpha_cond_pipe        "Heat Transfer Coefficient through condensed water film on pipe";
  SI.CoefficientOfHeatTransfer alpha_cond_plate       "Heat Transfer Coefficient through condensed water film on plates";

  SI.Length delta_pipe    "condensed water film thickness on pipe";
  SI.Length delta_plate   "condensed water film thickness on plates";

  //outlet state quantities
//  Modelica.Units.SI.Pressure p_out_water "pressure of medium exiting";
//  Modelica.Units.SI.SpecificEnthalpy h_out_water "enthaply of medium exiting";
//  Medium_Water.MassFraction Xi_out_water[Medium_Water.nXi] "mass fraction of medium exiting";

  ThermofluidStream.Interfaces.Inlet inlet_water(redeclare package Medium = Medium_Water) annotation (Placement(transformation(extent={{-120,-80},{-80,-40}}), iconTransformation(extent={{-120,-80},{-80,-40}})));
  ThermofluidStream.Interfaces.Outlet outlet_water(redeclare package Medium = Medium_Water) annotation (Placement(transformation(extent={{80,-80},{120,-40}}), iconTransformation(extent={{80,-80},{120,-40}})));

  ThermofluidStream.Topology.JunctionT1 junctionT1_1(redeclare package Medium = Medium_Water)                                        annotation (Placement(transformation(extent={{10,-70},{-10,-50}})));
  ThermofluidStream.Boundaries.Source source(redeclare package Medium = Medium_Water,
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

//protected
//  outer ThermofluidStream.DropOfCommons dropOfCommons;

initial equation
   dm_flow = 0.0;

equation

// Inflow and Wall resp. Tin Temperature
   T_gas_in = Medium.temperature(inlet.state);
   T_wall = heatPort.T;

// Enthalpy Inflow
   H_flow_in = h_in * m_flow;

// Density at inflow
   rho_in =Medium.density(inlet.state);

// Dynamic Viscosity at inflow
   eta_in =Medium.dynamicViscosity(inlet.state);

// Hollow Space Ratio
   eps = 1 - (Modelica.Constants.pi/4) * (d_out/H_ch);

// Mean Fluid Velocity in flow channel
   v_mean =inlet.m_flow/(rho_in*A_flow);

// Mean Fluid Velocity around pipe
   v_around = v_mean / eps;

// Reynolds Numbers
   Re_pipe = max(1.0,abs(v_around) * L_around * rho_in/eta_in);
   Re_rect = max(1.0,abs(v_mean)   * d_hyd    * rho_in/eta_in);

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
   lambda_F =Medium.thermalConductivity(inlet.state);

// Heat Capacity at inflow
   cp_F =Medium.specificHeatCapacityCp(inlet.state);

  alpha_air_pipe  = lambda_F/d_out * Nu_pipe;
  alpha_air_plate = lambda_F/d_hyd * Nu_plate;

  kA_pipe  = alpha_air_pipe  * A_HX_pipe_out;
  kA_plate = alpha_air_plate * A_HX_parallel;

  Q_flow_pipe  = kA_pipe * (T_wall - T_gas_in);
  Q_flow_plate = if noEvent(m_flow > 0) then cp_F * m_flow * (T_wall - T_gas_in) * (1 - exp( -kA_plate / (cp_F * m_flow))) else kA_plate*(T_wall - T_gas_in);

//------------------------------------------------------------------------
// Check if condensation can take place
//------------------------------------------------------------------------

  X_humid_in   = inlet.state.X[Medium.Water];
  X_vap_sat_in = Medium.Xsaturation(inlet.state);
  X_vap_in     = min(X_humid_in, X_vap_sat_in);
  X_liq_in     = X_humid_in - X_vap_in;

  m_flow_vap_in = X_vap_in * m_flow;

  p_vap_in = X_vap_in * R_vap / (R_da * (1-X_vap_in) + X_vap_in * R_vap) * p_in;

  p_vap_sat_wall = Medium.saturationPressureLiquid(T_wall);

  if noEvent(p_vap_in > p_vap_sat_wall + 1/1000) then

//------------------------------------------------------------------------
// Material Constants under current conditions
//------------------------------------------------------------------------

//  Gas Density
    rho_gas_wall = ThermofluidStream.Media.myMedia.Air.DryAirNasa.density_pT(p_in, T_wall);

//  Material Coefficients of liquid water at wall conditions
    rho_liq_wall = StandardWater.density_pT(p_in, T_wall);
    lambda_liq   = StandardWater.thermalConductivity(StandardWater.setState_pTX(p_in, T_wall));
    eta_liq      = StandardWater.dynamicViscosity(StandardWater.setState_pTX(p_in, T_wall));

//  Enthalpy of Condensation
    dh_vap = Medium.enthalpyOfVaporization(T_wall);

    cp_air_in =DryAirNasa.specificHeatCapacityCp(inlet.state);
    cp_vap_in =ThermofluidStream.Media.myMedia.IdealGases.SingleGases.H2O.specificHeatCapacityCp(inlet.state);
    cp_gas_in = (cp_air_in * (1 - X_vap_in) + cp_vap_in * X_vap_in) / (1-X_liq_in);

//  Material Transfer Coefficient
    beta_air_pipe  = alpha_air_pipe  / (cp_gas_in * rho_in) * Le_vap^(-2/3);  //TODO: wird nicht benötigt
    beta_air_plate = alpha_air_plate / (cp_gas_in * rho_in) * Le_vap^(-2/3);

//  X_vap_sat_wall = Medium.Xsaturation(Medium.setState_pTX(p_in, T_wall));

//------------------------------------------------------------------------
// Vapor Transfer to Condensation Area
//------------------------------------------------------------------------

//  Saturation Temperature at inflow Pressure
    T_sat = Medium.saturationTemperature(p_vap_in);

//  First Estimation of Condensation Mass Flow
//  ------------------------------------------
    T_cond_0 = (T_wall + T_sat) / 2.0;

//  Vapor Pressure on water film surface
    p_vap_sat_0 = Medium.saturationPressureLiquid(T_cond_0);  // TODO: Equal to p_vap_sat_wall ?

//  Condensation mass as diffusion flow to condensation area
    m_flow_cond_pipe_0  = alpha_air_pipe  / cp_gas_in * Le_vap^(-2/3) * A_HX_pipe_out * (p_vap_in - p_vap_sat_0) / (p_in - p_vap_sat_0);
    m_flow_cond_plate_0 = alpha_air_plate / cp_gas_in * Le_vap^(-2/3) * A_HX_parallel * (p_vap_in - p_vap_sat_0) / (p_in - p_vap_sat_0);

//  1st Iteration of Condensation Mass Flow //TODO
//  ---------------------------------------
    T_cond_pipe_1  = min(T_sat, max(T_wall, T_wall + (dh_vap / H_ch/ lambda_liq) * ((3.0 / (L_ch * sqrt(8.0)) * m_flow_cond_pipe_0)^4.0 * (eta_liq/(Modelica.Constants.g_n * rho_liq_wall * (rho_liq_wall-rho_gas_wall))))^(1.0/3.0)));
    T_cond_plate_1 = min(T_sat, max(T_wall, T_wall + (dh_vap / H_ch/ lambda_liq) * ((3.0 / (L_ch * sqrt(8.0)) * m_flow_cond_plate_0)^4.0 * (eta_liq/(Modelica.Constants.g_n * rho_liq_wall * (rho_liq_wall-rho_gas_wall))))^(1.0/3.0)));

//  Vapor Pressure at plates
    p_vap_sat_pipe_1  = Medium.saturationPressureLiquid(T_cond_pipe_1);
    p_vap_sat_plate_1 = Medium.saturationPressureLiquid(T_cond_plate_1);

//  Condensation mass as diffusion flow to condensation area
    m_flow_cond_pipe_1  = alpha_air_pipe  / cp_gas_in * Le_vap^(-2/3) * A_HX_pipe_out * (p_vap_in - p_vap_sat_pipe_1)  / (p_in - p_vap_sat_pipe_1);
    m_flow_cond_plate_1 = alpha_air_plate / cp_gas_in * Le_vap^(-2/3) * A_HX_parallel * (p_vap_in - p_vap_sat_plate_1) / (p_in - p_vap_sat_plate_1);

//  2nd Iteration of Condensation Mass Flow //TODO
//  ---------------------------------------
    if noEvent(T_cond_pipe_1 - T_cond_0 > 0.1) then
       T_cond_pipe  = min(T_sat, max(T_wall, T_wall + (dh_vap / H_ch/ lambda_liq) * ((3.0 / (L_ch * sqrt(8.0)) * m_flow_cond_pipe_1)^4.0 * (eta_liq/(Modelica.Constants.g_n * rho_liq_wall * (rho_liq_wall-rho_gas_wall))))^(1.0/3.0)));

       p_vap_sat_pipe  = Medium.saturationPressureLiquid(T_cond_pipe);

       m_flow_cond_pipe  = alpha_air_pipe  / cp_gas_in * Le_vap^(-2/3) * A_HX_pipe_out * (p_vap_in - p_vap_sat_pipe)  / (p_in - p_vap_sat_pipe);
    else
       T_cond_pipe      = T_cond_pipe_1;
       p_vap_sat_pipe   = p_vap_sat_pipe_1;
       m_flow_cond_pipe = m_flow_cond_pipe_1;
    end if;

    if noEvent(T_cond_plate_1 - T_cond_0 > 0.1) then
       T_cond_plate = min(T_sat, max(T_wall, T_wall + (dh_vap / H_ch/ lambda_liq) * ((3.0 / (L_ch * sqrt(8.0)) * m_flow_cond_plate_1)^4.0 * (eta_liq/(Modelica.Constants.g_n * rho_liq_wall * (rho_liq_wall-rho_gas_wall))))^(1.0/3.0)));

       p_vap_sat_plate = Medium.saturationPressureLiquid(T_cond_plate);

       m_flow_cond_plate = alpha_air_plate / cp_gas_in * Le_vap^(-2/3) * A_HX_parallel * (p_vap_in - p_vap_sat_plate) / (p_in - p_vap_sat_plate);
    else
       T_cond_plate      = T_cond_plate_1;
       p_vap_sat_plate   = p_vap_sat_plate_1;
       m_flow_cond_plate = m_flow_cond_plate_1;
    end if;

//  Total Condensation mass flow
    m_flow_cond = m_flow_cond_pipe + m_flow_cond_plate;

//  Ackermann Correction
    Phi_Ack_pipe  = m_flow_cond_pipe  * cp_gas_in / (A_HX_pipe_out * alpha_air_pipe);
    Phi_Ack_plate = m_flow_cond_plate * cp_gas_in / (A_HX_parallel * alpha_air_plate);

    zeta_Ack_pipe  = if noEvent(Phi_Ack_pipe  > 0.00001) then - Phi_Ack_pipe  / (exp(-Phi_Ack_pipe)  - 1) else 1.0;
    zeta_Ack_plate = if noEvent(Phi_Ack_plate > 0.00001) then - Phi_Ack_plate / (exp(-Phi_Ack_plate) - 1) else 1.0;

    alpha_air_pipe_Ack  = zeta_Ack_pipe  * alpha_air_pipe;
    alpha_air_plate_Ack = zeta_Ack_plate * alpha_air_plate;

//  Heat Flow when condensation
    Q_flow_cond_pipe  = m_flow_cond_pipe  * dh_vap + alpha_air_pipe_Ack  * A_HX_pipe_out * (T_gas_in - T_cond_pipe);
    Q_flow_cond_plate = m_flow_cond_plate * dh_vap + alpha_air_plate_Ack * A_HX_parallel * (T_gas_in - T_cond_plate);

    Q_flow_cond = Q_flow_cond_pipe + Q_flow_cond_plate;

    Q_flow = -Q_flow_cond;

//  Mean Heat Transfer Coefficient through condensed water film
    alpha_cond_pipe  = if noEvent(T_cond_pipe  - T_wall > 0.0001) then Q_flow_cond_pipe  / (A_HX_pipe_out * (T_cond_pipe  - T_wall)) else 1.0E+06;
    alpha_cond_plate = if noEvent(T_cond_plate - T_wall > 0.0001) then Q_flow_cond_plate / (A_HX_parallel * (T_cond_plate - T_wall)) else 1.0E+06;

//  Mean condensed water film thickness
    delta_pipe  = A_HX_pipe_out * lambda_liq * (T_cond_pipe  - T_wall)/Q_flow_cond_pipe;
    delta_plate = A_HX_parallel * lambda_liq * (T_cond_plate - T_wall)/Q_flow_cond_plate;

//  Specific condensed water enthalpy
    h_water_pipe  = StandardWater.specificEnthalpy_pT(p_in, (T_wall+T_cond_pipe) /2);
    h_water_plate = StandardWater.specificEnthalpy_pT(p_in, (T_wall+T_cond_plate)/2);

//  Enthalpy Flow of water
    H_flow_water_pipe  = h_water_pipe  * m_flow_cond_pipe;
    H_flow_water_plate = h_water_plate * m_flow_cond_plate;

    H_flow_water = H_flow_water_pipe + H_flow_water_plate;

    T_water = if noEvent(m_flow_cond > 0) then ((T_cond_pipe + T_wall)/2 * m_flow_cond_pipe + (T_cond_plate + T_wall)/2 * m_flow_cond_plate) / m_flow_cond else T_wall;

//  Mean Heat Transfer Coefficient through the water film on vertical wall
//    alpha_liq_vert_0 = 2/3*sqrt(2) * (dh_vap * Modelica.Constants.g_n * rho_liq_wall * (rho_liq_wall-rho_gas_wall) * lambda_liq^3 /
//                                     (eta_liq * (T_cond_0 - T_wall) * H_ch)) ^(1/4);

//  1st Iteration of Condensation Mass Flow
//  ---------------------------------------

//  Ackermann Correction
/*    Phi_Ack  = Le_vap^(-2/3) * (p_vap_in - p_vap_sat_0) / (p_in - p_vap_sat_0);
    zeta_Ack = if noEvent(Phi_Ack > 0.00001) then - Phi_Ack / (exp(-Phi_Ack) - 1) else 1.0;

    alpha_air_pipe_Ack  = zeta_Ack * alpha_air_pipe;
    alpha_air_plate_Ack = zeta_Ack * alpha_air_plate;*/

//  Condensation Temperature on plate from heat flow balance
//    T_cond_plate_1 = min(T_sat, max(T_wall, (dh_vap * m_flow_cond_plate_0 / A_HX_parallel + alpha_air_plate_Ack * T_gas_in + alpha_liq_vert_0 * T_wall) /
//                                            (alpha_air_plate_Ack + alpha_liq_vert_0)));

//  Vapor Pressure at plates
//    p_vap_sat_plate_1 = Medium.saturationPressureLiquid(T_cond_plate_1);

//  Mean Heat Transfer Coefficient through the water film on vertical wall
//    alpha_liq_vert_1 = 2/3*sqrt(2) * (dh_vap * Modelica.Constants.g_n * rho_liq_wall * (rho_liq_wall-rho_gas_wall) * lambda_liq^3 /
//                                     (eta_liq * (T_cond_plate_1 - T_wall) * H_ch)) ^(1/4);

//  Condensation mass as diffusion flow to condensation area
//    m_flow_cond_plate_1 = max(0.0, alpha_air_plate_Ack / cp_gas_in * Le_vap^(-2/3) * A_HX_parallel * (p_vap_in - p_vap_sat_plate_1) / (p_in - p_vap_sat_plate_1));

    // Theoretic Mean Film Thickness
//    delta_film_mean_1 = ((3 * eta_liq * m_flow_cond_0) / (Modelica.Constants.g_n * rho_liq_wall * (rho_liq_wall-rho_gas_wall) * L_cond))^(1/3);

//    Q_flow_plate_cond = alpha_liq_vert_1 * A_HX_parallel * (T_wall - T_cond_plate_1);

//    m_flow_vap_out = m_flow_vap_in - m_flow_cond_plate_1;

// ------------------------------------------------------------------------
   else    //  Heat Flow of Humid Air without Condensation
// ------------------------------------------------------------------------

    rho_gas_wall = 0.0;
    rho_liq_wall = 0.0;
    lambda_liq   = 0.0;
    eta_liq      = 0.0;

    dh_vap = 0.0;

    cp_air_in = 0.0;
    cp_gas_in = 0.0;
    cp_vap_in = 0.0;

    beta_air_pipe  = 0.0;
    beta_air_plate = 0.0;

    T_sat          = 0.0;
    T_cond_0       = 0.0;
    T_cond_pipe_1  = 0.0;
    T_cond_plate_1 = 0.0;
    T_cond_pipe    = 0.0;
    T_cond_plate   = 0.0;

    p_vap_sat_0       = 0;
    p_vap_sat_pipe_1  = 0;
    p_vap_sat_plate_1 = 0;
    p_vap_sat_pipe    = 0;
    p_vap_sat_plate   = 0;
//    p_vap_sat_2 = 0;

    m_flow_cond_pipe_0 = 0.0;
    m_flow_cond_pipe_1 = 0.0;
    m_flow_cond_pipe   = 0.0;
    m_flow_cond        = 0.0;
    //    m_flow_cond_1 = 0;
//    m_flow_cond_2 = 0;
//    m_flow_cond_cell = 0;

    m_flow_cond_plate_0 = 0.0;
    m_flow_cond_plate_1 = 0.0;
    m_flow_cond_plate   = 0.0;

//    m_flow_vap_out   = 0.0;

    Phi_Ack_pipe  = 0.0;
    Phi_Ack_plate = 0.0;
    zeta_Ack_pipe  = 0.0;
    zeta_Ack_plate = 0.0;
    alpha_air_pipe_Ack  = 0.0;
    alpha_air_plate_Ack = 0.0;
//    alpha_liq_vert_0 = 0.0;
//    alpha_liq_vert_1 = 0.0;

    alpha_cond_pipe  = 0.0;
    alpha_cond_plate = 0.0;

    delta_pipe  = 0.0;
    delta_plate = 0.0;

    Q_flow_cond_pipe  = 0.0;
    Q_flow_cond_plate = 0.0;
    Q_flow_cond       = 0.0;

    Q_flow = Q_flow_pipe + Q_flow_plate;

    h_water_pipe  = 0.0;
    h_water_plate = 0.0;

    H_flow_water_pipe  = 0.0;
    H_flow_water_plate = 0.0;
    H_flow_water       = 0.0;

    T_water = 290.0;

  end if;

//  outlet_liquid_water.

//------------------------------------------------------------------------
// Final Air and Water Mass Flow
//------------------------------------------------------------------------

   der(dm_flow) = (-m_flow_cond - dm_flow) / dt;

   m_flow_out = m_flow + dm_flow;

   m_flow_humid_out = X_humid_in * m_flow + dm_flow;   // TODO m_flow = m_flow_in

//------------------------------------------------------------------------
// Final Outflow Humdity
//------------------------------------------------------------------------

   X_humid_out = if noEvent(m_flow_out > 0) then m_flow_humid_out / m_flow_out else X_humid_in;

   Xi_out = {X_humid_out};

//------------------------------------------------------------------------
// Final Outflow Heat
//------------------------------------------------------------------------

// Outflow Specific Enthalpy
// -------------------------
//   Q_flow_cond = Q_flow_cond_pipe + Q_flow_cond_plate;

//   H_flow_water = H_flow_water_pipe + H_flow_water_plate;

   H_flow_out = H_flow_in + Q_flow - H_flow_water;

   h_out = if noEvent(m_flow_out > 0 and H_flow_out >0) then H_flow_out / m_flow_out else h_in;

//   h_out_water = if noEvent(m_flow_out > 0) then -H_flow_water / dm_flow else 0.1;

   //dm_flow     = -m_flow_cond_plate_1; // bringt Fehler

   //outlet_water.r = 0.0;

// -------------------------------------------------------------------------------
// Pressure Drop against a crosswise pipe
//--------------------------------------------------------------------------------

// Dynamic Fluid Pressure around pipe
   p_dyn = rho/2 * v_around^2;

// Pressure Drop from drag force of a cylinder
   dp = - cd_cyl * d_out/H_ch * p_dyn;

  connect(junctionT1_1.outlet, outlet_water) annotation (Line(
      points={{10,-60},{100,-60}},
      color={28,108,200},
      thickness=0.5));
  connect(inlet_water, junctionT1_1.inletB) annotation (Line(
      points={{-100,-60},{-40,-60},{-40,-88},{0,-88},{0,-70}},
      color={28,108,200},
      thickness=0.5));
  connect(Temperature_Water_condensed.y, source.T0_var) annotation (Line(points={{13,40},{0,40},{0,22}}, color={0,0,127}));
  connect(source.outlet, condensedWaterFlowSensor.inlet) annotation (Line(
      points={{0,10},{0,-10}},
      color={28,108,200},
      thickness=0.5));
  connect(condensedWaterFlowSensor.outlet, junctionT1_1.inletA) annotation (Line(
      points={{0,-30},{0,-50}},
      color={28,108,200},
      thickness=0.5));
  connect(MassFlow_condensed.y, condensedWater_PressureController.u_s) annotation (Line(points={{-69,40},{-61,40}}, color={0,0,127}));
  connect(condensedWater_PressureController.u_m, condensedWaterFlowSensor.value_out) annotation (Line(points={{-50,29},{-50,-40},{-6,-40},{-6,-30}}, color={0,0,127}));
  connect(condensedWater_PressureController.y, source.p0_var) annotation (Line(points={{-39,40},{-6,40},{-6,22}}, color={0,0,127}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
        Text(
          extent={{-68,-72},{70,-100}},
          textColor={0,0,0},
          textString="Moist Air (condensing)")}),
    Diagram(coordinateSystem(preserveAspectRatio=false)),
    Documentation(info="<html>
    For a single conduction element this component models
    <ul>
    <li>the temperature drop or rise of the fluid along the cross channel</li>
    <li>the heat flow between fluid and pipe wall resp. fins</li>
    <li>the pressure drop along the cross channel.</li>
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
    The Nusselt number for flows around pipes is calculated according to Gnielinksi.<br>
    as documented with formula (3.211) in the text book from H.D. Baehr and K. Stephan: \"W&auml;rme- und Stoff&uuml;bertragung\", 7th edition.<br>
    <p>
    The total heat exchange rate is the sum of heat exchange to the pipe and to the fins.<br>
    The heat exchange rate to the fins is calulated at the beginning of the fins. Along the length of the fins an exponential convergence to the wall temperature is assumed.
    <p>
    <em>Condensation</em>
    <p>
    Condensation can take place when the vapor pressure of the air exceeds the saturation vapor pressure at wall temperature.<br>
    The condensation mass flow is mainly determined by the transport of the distributed vapor to the condensation areas.
    This is a diffusion process which is driven by the vapor pressure difference as expressed in formula (4.31)<br>
    The condensed water forms a water film on the surfaces, which affects the heat exchange rate as described in formula (4.30).<br>
    All the effects are described in chapter 4.14 \Einfluss nicht kondensierbarer Gase\" in the named text book.
    
</html>"));
end ConductionElement_Cross_CondensingAir;
