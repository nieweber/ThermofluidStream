within ThermofluidStream.HeatExchangersPhysical.ConductionElements;
model ConductionElement_Pipe_2Phase "ConductionElement for one two-phase fluid in pipe"
  extends PipeIcon;
  extends Internal.PartialConductionCell(
                                redeclare package Medium =
        ThermofluidStream.Media.myMedia.Interfaces.PartialTwoPhaseMedium);  // Refrigerant

  import Modelica.Units.SI;

  parameter SI.Length d_hyd = 1   "Hydraulic Diameter";
  parameter SI.Length L_flow = 1  "Length of flow";

  final parameter SI.Area A_flow = Modelica.Constants.pi * (d_hyd/2)^2;
  final parameter SI.Area A_inside =  Modelica.Constants.pi * d_hyd * L_flow "inner lateral area";

  constant SI.Duration dt_cond = 0.1  "time constant of condensation";

  constant SI.SpecificEnthalpy h_tol = 1.0 "Tolerance to separate regions";
  constant SI.SpecificEnthalpy h_tol_evap = 100 "Tolerance to separate regions";

  constant SI.ReynoldsNumber Re_crit_low =  1000  "lower limit of transition phase laminar to turbulent";
  constant SI.ReynoldsNumber Re_crit_high = 2000  "upper limit of transition phase laminar to turbulent";

  SI.Density rho_in      "2-phase density at inlet";
  SI.Density rho_bubble  "density in liquid phase at inlet pressure";
  SI.Density rho_dew     "density at dew point at inlet pressure";

  Real phase_state;

  Medium.SaturationProperties sat_in;

  SI.Temperature T_in    "inflow temperature";
  SI.Temperature T_wall  "wall temperature";
  SI.Temperature T_sat   "condensation temperature";

  SI.SpecificEnthalpy dh_vap  "enthalpy of vaporization";

  SI.Velocity v_mean     "mean fluid velocity";

  SI.DynamicViscosity mu      "dynamic viscosity";
  SI.DynamicViscosity mu_bubble  "dynamic viscosity at bubble point";
  SI.DynamicViscosity mu_dew     "dynamic viscosity at dew point";

  SI.ReynoldsNumber Re           "Reynolds number";
  SI.ReynoldsNumber Re_liq       "Reynolds number in liquid phase";
  SI.ReynoldsNumber Re_gas       "Reynolds number in gas phase";
  SI.ReynoldsNumber Re_liq_part  "Reynolds number of liquid partition";
  SI.ReynoldsNumber Re_gas_part  "Reynolds number of gas partition";

  SI.PrandtlNumber Pr_F  "Prandtl Number at inlet temperature";
  SI.PrandtlNumber Pr_W  "Prandtl Number at wall temperature";

  SI.PrandtlNumber Pr_liq  "Prandtl Number in liquid phase";
  SI.PrandtlNumber Pr_gas  "Prandtl Number in gas phase";

  SI.PecletNumber Pe      "Peclet number";
  SI.PecletNumber Pe_liq  "Peclet number in liquid phase";
  SI.PecletNumber Pe_gas  "Peclet number in gas phase";

  Real L_norm       "Normalized pipe length";
  Real L_norm_liq   "Normalized pipe length in liquid phase";
  Real L_norm_gas   "Normalized pipe length in gas phase";

  SI.NusseltNumber Nu       "Nusselt Number";
  SI.NusseltNumber Nu_liq   "Nusselt Number in liquid phase";
  SI.NusseltNumber Nu_gas   "Nusselt Number in gas phase";

  Real zeta       "Pipe friction zeta value";
  Real zeta_liq   "Pipe friction zeta value in liquid phase";
  Real zeta_gas   "Pipe friction zeta value in gas phase";

  Real pipeFric_liq       "Pipe friction lambda value in liquid phase";
  Real pipeFric_gas       "Pipe friction lambda value in gas phase";

  Real pipeFric_liq_lam   "Pipe friction lambda value in laminar liquid phase";
  Real pipeFric_gas_lam   "Pipe friction lambda value in laminar gas phase";

  Real pipeFric_liq_turb  "Pipe friction lambda value in turbulent liquid phase";
  Real pipeFric_gas_turb  "Pipe friction lambda value in turbulent gas phase";

  Real K_fac  "K factor for liquid fluids";

  SI.CoefficientOfHeatTransfer alpha_mean  "Mean Coefficient of Heat transfer";
  SI.CoefficientOfHeatTransfer alpha_gas   "Heat Transfer between gas and condensed film or wall";
  SI.CoefficientOfHeatTransfer alpha_cond  "Heat Transfer Coefficient by condensation";
  SI.CoefficientOfHeatTransfer alpha_liq   "Heat Transfer Coefficient in liquid phase";

  SI.ThermalConductance kA       "Thermal conductance fluid -> wall";
  SI.ThermalConductance kA_cond  "Thermal conductance of condensate";

  SI.ThermalConductivity lambda_F       "Thermal Conductivity of Fluid";
  SI.ThermalConductivity lambda_bubble  "Thermal Conductivity at bubble point";
  SI.ThermalConductivity lambda_dew     "Thermal Conductivity at dew point";

  SI.SpecificHeatCapacityAtConstantPressure cp_bubble  "specific heat capacity at bubble point";
  SI.SpecificHeatCapacityAtConstantPressure cp_dew     "specific heat capacity at dew point";

  Real X_gas     "gas partition in cell";
  Real X_liq     "liquid partition in cell";

  Real X_gas_in  "gas partition at inlet";
  Real X_liq_in  "liquid partition at inlet";

  Real X_gas_inflow  "gas partition of cell mass inflow";
  Real X_liq_inflow  "liquid partition of cell mass inflow";

  Real X_liq_lam   "laminar portion in liquid phase";
  Real X_liq_turb  "turbulent portion in liquid phase";

  Real X_gas_lam   "laminar portion in gas phase";
  Real X_gas_turb  "turbulent portion in gas phase";

  Real C      "flow type dependent constant of Lockhart-Martinelli parameter";
  Real X2     "Lockhart-Martinelli parameter";
  Real Phi_2  "weighting function of Lockhart-Martinelli parameter";

  SI.Area A_liq   "inner lateral area covered with liquid";
  SI.Area A_cond  "inner lateral area ready for condensation";

  SI.Length delta_film_mean     "Mean Thickness of condensed film in pipe";

  SI.Mass m_gas  "gas mass in cell";

  SI.MassFlowRate m_flow_gas_in       "mass flow rate of gas partition at inlet";
  SI.MassFlowRate m_flow_cond         "mass flow rate of flowing condensate";
  SI.MassFlowRate m_flow_cond_avail   "available condensation rate due to existing gas and condensation rate";
  SI.MassFlowRate m_flow_cond_pot     "potential mass flow rate of flowing condensate (Nusselt Formula)";
  SI.MassFlowRate m_flow_cond_equlib  "mass flow rate of flowing condensate in equilibrium (Nusselt Formula)";

  Real m_flow_density(unit="kg/(s.m2)")      "mass flow density";
  Real m_flow_density_liq(unit="kg/(s.m2)")  "mass flow density of liquid partition";
  Real m_flow_density_gas(unit="kg/(s.m2)")  "mass flow density of gas partition";

  SI.HeatFlowRate Q_flow2           "heat flow";
  SI.HeatFlowRate Q_flow_cond       "heat flow from condensate film to wall";
  SI.HeatFlowRate Q_flow_cond_pot   "potential heat flow from condensate film to wall";
  SI.HeatFlowRate Q_flow_gas_film   "heat flow from gas to condensate film";
  SI.HeatFlowRate Q_flow_gas_film2  "heat flow from gas to condensate film";
  SI.HeatFlowRate Q_flow_gas_wall   "heat flow from gas to wall";
  SI.HeatFlowRate Q_flow_gas_wall2  "heat flow from gas to wall";
  SI.HeatFlowRate Q_flow_liq_wall   "heat flow from liquid to wall";

  SI.SpecificEnthalpy h_bubble  "specific enthalpy at bubble point";
  SI.SpecificEnthalpy h_dew     "specific enthalpy at dew point";

  Real dpdz_liq(unit="Pa/m")    "pressure gradient of liquid phase";
  Real dpdz_gas(unit="Pa/m")    "pressure gradient of gas phase";
  Real dpdz_fric(unit="Pa/m")   "pressure gradient by friction";

initial equation
//  kA = 1;
//  v_mean = 1;

equation

  T_in   = Medium.temperature(inlet.state);
  T_wall = heatPort.T;

  sat_in = Medium.setSat_p(p_in);   //Refrigerant

  h_dew    = Medium.dewEnthalpy(sat_in);      //Refrigerant
  h_bubble = Medium.bubbleEnthalpy(sat_in);   //Refrigerant
  dh_vap   = h_dew - h_bubble;

  T_sat = Medium.saturationTemperature(p_in);

  rho_in     = Medium.density(inlet.state);
  rho_bubble = Medium.bubbleDensity(sat_in);
  rho_dew    = Medium.dewDensity(sat_in);

  lambda_F      = Medium.thermalConductivity(inlet.state);
  lambda_bubble = Medium.thermalConductivity(Medium.setState_pTX(p_in, T_sat - 1, Xi_in));
  lambda_dew    = Medium.thermalConductivity(Medium.setState_pTX(p_in, T_sat + 1, Xi_in));

  cp_bubble = Medium.specificHeatCapacityCp(Medium.setState_pTX(p_in, T_sat - 1, Xi_in));
  cp_dew    = Medium.specificHeatCapacityCp(Medium.setState_pTX(p_in, T_sat + 1, Xi_in));

  mu        = Medium.dynamicViscosity(inlet.state);
  mu_bubble = Medium.dynamicViscosity(Medium.setState_pTX(p_in, T_sat - 1, Xi_in));
  mu_dew    = Medium.dynamicViscosity(Medium.setState_pTX(p_in, T_sat + 1, Xi_in));

  X_gas_in = max(0, min(1, (h_in - h_bubble)/(h_dew - h_bubble)));
  X_liq_in = 1 - X_gas_in;

  A_cond = A_inside * X_gas_in;
  A_liq  = A_inside * X_liq_in;

  X_gas_inflow = X_gas_in;
  X_liq_inflow = X_liq_in;

  m_flow_gas_in = m_flow * X_gas_in;

  X_gas = max(0, min(1, (h - h_bubble)/(h_dew - h_bubble)));
  X_liq = 1 - X_gas;
  m_gas = M * X_gas;

  v_mean = inlet.m_flow/(rho_in*A_flow);

  Re = max(1.0, abs(v_mean) * d_hyd * rho_in/mu);

//-----------------------------------------------------------------------------------------
  if noEvent(T_wall >= T_sat and h_in <= h_dew) then  // heating the pipe and evaporation
//-----------------------------------------------------------------------------------------

    L_norm = 0;
    kA_cond = 0;
    m_flow_cond = 0;
    m_flow_cond_pot = 0;
    m_flow_cond_avail = 0;
    m_flow_cond_equlib = 0;
    delta_film_mean = 0;
    alpha_cond = 0;
    Q_flow_cond = 0;
    Q_flow_cond_pot = 0;
    Q_flow_gas_film = 0;
    Q_flow_gas_film2 = 0;
    Q_flow_gas_wall = 0;
    Q_flow_gas_wall2 = 0;

    Pr_F = 0;
    Pe = 0;
    Q_flow_liq_wall = 0;
//    dp = 0;

    zeta = 0;

/*  Heat Transfer Coefficient in Liquid Phase */

    if noEvent(h_in < h_bubble) then  // inflow in liquid phase

      phase_state = 11;

      Pr_liq = Medium.prandtlNumber(inlet.state);

      if noEvent(T_wall > T_sat + 1) then
        Pr_W  = Medium.prandtlNumber(Medium.setState_pTX(p_in, T_wall, Xi_in));
      else
        Pr_W  = Medium.prandtlNumber(Medium.setState_pTX(p_in, T_sat + 1, Xi_in));
      end if;

      K_fac = (Pr_liq/Pr_W)^0.11;

      Re_liq = Re;
      Re_gas = Re;

    else  // inflow in 2-phase flow

      phase_state = 12;

      Pr_liq = Medium.prandtlNumber(Medium.setState_pTX(p_in, T_sat - 1, Xi_in));

      Pr_W = 0;

      K_fac = 1;

      Re_liq = max(1.0,abs(v_mean) * d_hyd * rho_bubble/mu_bubble);
      Re_gas = max(1.0,abs(v_mean) * d_hyd * rho_dew/mu_dew);

    end if;

    Pr_gas = Medium.prandtlNumber(Medium.setState_pTX(p_in, T_sat + 1, Xi_in));

    Pe_liq = Re_liq * Pr_liq;
    Pe_gas = Re_gas * Pr_gas;

    L_norm_liq = L_flow / (d_hyd * Pe_liq);
    L_norm_gas = L_flow / (d_hyd * Pe_gas);

    zeta_liq = (0.78 * log(Re_liq) - 1.5)^(-2);
    zeta_gas = (0.78 * log(Re_gas) - 1.5)^(-2);

    if noEvent(Re_liq < 2300) then  // Laminar Flow

      Nu_liq = 3.657 / tanh(2.264*L_norm_liq^(1./3.) + 1.7*L_norm_liq^(2./3.)) + 0.0499/L_norm_liq*tanh(L_norm_liq);    // Baehr-Stephan (3.250)
      //Nu_liq = (49.37 + (1.615*(Pe_liq * d_hyd/L_flow)^(1/3) -0.7)^3)^(1/3) * K_fac;                  // schweizer-fn

    else  // Turbulent Flow

      Nu_liq = (zeta_liq/8)*Re_liq*Pr_liq / (1+12.7*(Pr_liq^(2/3)-1)*sqrt(zeta_liq/8)) * (1+(d_hyd/L_flow)^(2/3)) * K_fac;           // Baehr-Stephan (3.260)
      // Nu_liq = zeta_liq/8*(Re_liq-1000)*Pr_liq / (1+12.7*(Pr_liq^(2/3)-1)*sqrt(zeta_liq/8)) * (1+(d_hyd/L_flow)^(2/3)) * K_fac;       // schweizer-fn

    end if;

    if noEvent(Re_gas < 2300) then  // Laminar Flow

      Nu_gas = 3.657 / tanh(2.264*L_norm_gas^(1./3.) + 1.7*L_norm_gas^(2./3.)) + 0.0499/L_norm_gas*tanh(L_norm_gas);    // Baehr-Stephan (3.250)
      //Nu_gas = (49.37 + (1.615*(Pe_gas * d_hyd/L_flow)^(1/3) -0.7)^3)^(1/3) * K_fac;                  // schweizer-fn

    else  // Turbulent Flow

      Nu_gas = (zeta_gas/8)*Re_gas*Pr_gas / (1+12.7*(Pr_gas^(2/3)-1)*sqrt(zeta_gas/8)) * (1+(d_hyd/L_flow)^(2/3)) * K_fac;           // Baehr-Stephan (3.260)
      // Nu_gas = zeta_gas/8*(Re_gas-1000)*Pr_gas / (1+12.7*(Pr_gas^(2/3)-1)*sqrt(zeta_liq/8)) * (1+(d_hyd/L_flow)^(2/3)) * K_fac;       // schweizer-fn

    end if;

    Nu = 0;

    if noEvent(h_in <= h_bubble) then  // inflow in liquid phase

      alpha_liq = lambda_F/d_hyd * Nu_liq;

    else  // inflow in 2-phase flow

      alpha_liq = lambda_bubble/d_hyd * Nu_liq;

    end if;

    alpha_gas = lambda_dew/d_hyd * Nu_gas;

    // Baehr-Stephan (4.160)
    alpha_mean = alpha_liq * (X_liq_inflow^0.01 * (X_liq_inflow + 1.2*X_gas_inflow^0.4 * (rho_bubble/rho_dew)^0.37)^(-2.2) +
                              X_gas_inflow^0.01 * (alpha_gas/alpha_liq * (1 + 8*X_liq_inflow^0.7 * (rho_bubble/rho_dew)^0.67))^(-2))^(-0.5);

    kA = alpha_mean * A_inside;

    if noEvent(h_in <= h_bubble) then  // inflow in liquid phase

      Q_flow  = if noEvent(m_flow > 0.001) then cp_bubble * m_flow * (T_wall - T_in) * (1 - exp( -kA / (cp_bubble * m_flow))) else kA*(T_wall - T_in);
      Q_flow2 = kA*(T_wall - T_in);

    else  // inflow in 2-phase flow

      Q_flow  = min(m_flow*(h_dew-h_in+h_tol_evap), kA*(T_wall - T_in));
      Q_flow2 = if noEvent(m_flow > 0.001) then cp_bubble * m_flow * (T_wall - T_in) * (1 - exp( -kA / (cp_bubble * m_flow))) else kA*(T_wall - T_in);

    end if;

//---------------------------------------------------------------------------------------------------------
elseif noEvent(T_wall <= T_sat and h_in > h_bubble + h_tol) then   // cooling the pipe with (potential) condensation
//---------------------------------------------------------------------------------------------------------

    phase_state = if noEvent(h_in > h_dew) then 23 else 22;

      m_flow_cond_avail = m_gas / dt_cond;

      m_flow_cond_equlib = 2.28714 * L_flow * ((d_hyd * (T_sat - T_wall) * lambda_bubble) / dh_vap)^(3/4) * ((Modelica.Constants.g_n * rho_bubble * (rho_bubble-rho_dew)) / mu_bubble)^(1/4);

      m_flow_cond_pot = min(m_flow_cond_avail, m_flow_cond_equlib);

      Q_flow_cond_pot = - dh_vap * m_flow_cond_pot;

    if noEvent(h_in <= h_dew) then
      Pr_F = Medium.prandtlNumber(Medium.setState_pTX(p_in, T_sat + 1, Xi_in));
    else
      Pr_F = Medium.prandtlNumber(inlet.state);
    end if;

    Pe   = Re * Pr_F;
    K_fac = 1;

    zeta = (0.79 * log(Re) - 1.64)^(-2);

    if noEvent(Re > 2300) then // turbulent
      Nu = zeta/8*(Re-1000)*Pr_F / (1+12.7*(Pr_F^(2/3)-1)*sqrt(zeta/8)) * (1+(d_hyd/L_flow)^(2/3)) * K_fac;
    else // laminar
      Nu = (49.37 + (1.615*(Pe*d_hyd/L_flow)^(1/3) -0.7)^3)^(1/3) * K_fac;
    end if;

    alpha_gas = lambda_F/d_hyd * Nu;

    Q_flow_gas_film2 = alpha_gas * (T_sat - T_in) * A_inside;
    Q_flow_gas_film = if noEvent(m_flow > 0.001) then cp_dew * m_flow * (T_sat - T_in) * (1 - exp( -(alpha_gas * A_inside) / (cp_dew * m_flow))) else alpha_gas*A_inside * (T_sat - T_in);

    Q_flow_gas_wall2 = alpha_gas * (T_wall - T_in) * A_inside;
    Q_flow_gas_wall = if noEvent(m_flow > 0.001) then cp_dew * m_flow * (T_wall - T_in) * (1 - exp( -(alpha_gas * A_inside) / (cp_dew * m_flow))) else alpha_gas*A_inside * (T_wall - T_in);

    if noEvent(Q_flow_gas_film <= Q_flow_cond_pot) then // no condensation
      Q_flow_cond = 0;
      Q_flow = Q_flow_gas_wall;
      m_flow_cond = 0;
      kA_cond = 0;
      alpha_cond = 0;
      delta_film_mean = 0;

    else  // condensation
      Q_flow_cond = Q_flow_cond_pot - Q_flow_gas_film;
      m_flow_cond = -Q_flow_cond / dh_vap;
      kA_cond = -Q_flow_cond / (T_sat - T_wall);
      alpha_cond = kA_cond / A_inside;
      delta_film_mean = lambda_bubble / alpha_cond;

      Q_flow = Q_flow_cond + Q_flow_gas_film;
    end if;

    kA = -Q_flow / max(1.E-5, abs(T_in - T_wall));
    alpha_mean = kA / A_inside;

    Pr_W = 0;
    Pr_liq = 0;
    Pr_gas = 0;
    Re_liq = 0;
    Re_gas = 0;
    Pe_liq = 0;
    Pe_gas = 0;
    L_norm = 0;
    L_norm_liq = 0;
    L_norm_gas = 0;
    zeta_liq = 0;
    zeta_gas = 0;
    Nu_liq = 0;
    Nu_gas = 0;
    alpha_liq = 0;
    Q_flow2 = 0;
    Q_flow_liq_wall = 0;
//    dp = 0;

// -------------------------------------------------------------------------------
else  // pure liquid or pure gas:  (h_in <= h_bubble + h_tol) or (h_in > h_dew)
//--------------------------------------------------------------------------------

    Pr_liq = 0;
    Pr_gas = 0;
    Re_liq = 0;
    Re_gas = 0;
    Pe_liq = 0;
    Pe_gas = 0;
    L_norm_liq = 0;
    L_norm_gas = 0;
    zeta_liq = 0;
    zeta_gas = 0;
    Nu_liq = 0;
    Nu_gas = 0;
    kA_cond = 0;
    m_flow_cond = 0;
    m_flow_cond_pot = 0;
    m_flow_cond_avail = 0;
    m_flow_cond_equlib = 0;
    delta_film_mean = 0;
    alpha_liq = 0;
    alpha_gas = 0;
    alpha_cond = 0;
    Q_flow_cond = 0;
    Q_flow_cond_pot = 0;
    Q_flow_liq_wall = 0;
    Q_flow_gas_film = 0;
    Q_flow_gas_film2 = 0;
    Q_flow_gas_wall = 0;
    Q_flow_gas_wall2 = 0;

    phase_state = if noEvent(h_in < (h_bubble + h_dew)/2) then 21 else 13;

    Pr_F  = Medium.prandtlNumber(inlet.state);

    if noEvent(T_wall < T_sat - 1) then
      Pr_W  = Medium.prandtlNumber(Medium.setState_pTX(p_in, T_wall, Xi_in));
      K_fac = (Pr_F/Pr_W)^0.11;
    elseif noEvent(T_wall > T_sat + 1) then
      Pr_W  = Medium.prandtlNumber(Medium.setState_pTX(p_in, T_wall, Xi_in));
      K_fac =(Pr_F/Pr_W)^0.11;
    else
      Pr_W  = 0;
      K_fac = 1;
    end if;

    Pe = Re * Pr_F;

    L_norm = L_flow / (d_hyd * Pe);

    zeta = (0.78 * log(Re) - 1.5)^(-2);       // Baehr-Stephan
    //zeta = (0.79 * log(Re) - 1.64)^(-2);    // schweizer-fn

    if noEvent(Re < 2300) then  // Laminar Flow

      Nu = 3.657 / tanh(2.264*L_norm^(1./3.) + 1.7*L_norm^(2./3.)) + 0.0499/L_norm*tanh(L_norm);    // Baehr-Stephan
      //Nu = (49.37 + (1.615*(Pe*d_hyd/L_flow)^(1/3) -0.7)^3)^(1/3) * K_fac;                        // schweizer-fn

    else  // Turbulent Flow

      Nu = (zeta/8)*Re*Pr_F / (1+12.7*(Pr_F^(2/3)-1)*sqrt(zeta/8)) * (1+(d_hyd/L_flow)^(2/3)) * K_fac;           // Baehr-Stephan
      // Nu = zeta/8*(Re-1000)*Pr_F / (1+12.7*(Pr_F^(2/3)-1)*sqrt(zeta/8)) * (1+(d_hyd/L_flow)^(2/3)) * K_fac;       // schweizer-fn

    end if;

    alpha_mean = lambda_F/d_hyd * Nu;

    kA = alpha_mean * A_inside;

    if noEvent(h_in < (h_bubble + h_dew)/2) then  // liquid flow
      Q_flow  = if noEvent(m_flow > 0.001) then cp_bubble * m_flow * (T_wall - T_in) * (1 - exp( -kA / (cp_bubble * m_flow))) else kA*(T_wall - T_in);
    else
      Q_flow  = if noEvent(m_flow > 0.001) then cp_dew * m_flow * (T_wall - T_in) * (1 - exp( -kA / (cp_dew * m_flow))) else kA*(T_wall - T_in);
    end if;

    Q_flow2 = kA*(T_wall - T_in);

//    dp = 0;

  end if;

// -------------------------------------------------------------------------------
// Pressure Drop
//--------------------------------------------------------------------------------

// Mass Flow Density
// -----------------
   m_flow_density = m_flow/A_flow;
   m_flow_density_liq = m_flow_density * X_liq;
   m_flow_density_gas = m_flow_density * X_gas;

// Reynolds Numbers of Flow Partitions
// -----------------------------------
   Re_liq_part = max(1.0, abs(m_flow_density_liq) * d_hyd / mu_bubble);
   Re_gas_part = max(1.0, abs(m_flow_density_gas) * d_hyd / mu_dew);

// Laminar/turbulent partitions
// ----------------------------
   (X_liq_lam, X_liq_turb) =Internal.portion(Re_crit_low, Re_crit_high, Re_liq_part);
   (X_gas_lam, X_gas_turb) =Internal.portion(Re_crit_low, Re_crit_high, Re_gas_part);

// Friction Factors
// ----------------
   pipeFric_liq_lam  = 64/Re_liq_part;
   pipeFric_gas_lam  = 64/Re_gas_part;

   pipeFric_liq_turb = 0.3164/(Re_liq_part^0.25);
   pipeFric_gas_turb = 0.3164/(Re_gas_part^0.25);

   pipeFric_liq = X_liq_lam * pipeFric_liq_lam + X_liq_turb * pipeFric_liq_turb;
   pipeFric_gas = X_gas_lam * pipeFric_gas_lam + X_gas_turb * pipeFric_gas_turb;

   //pipeFric_liq = transitionFunction(pipeFric_liq_lam, pipeFric_liq_turb, Re_crit_low, Re_crit_high, Re_liq_part);
   //pipeFric_gas = transitionFunction(pipeFric_gas_lam, pipeFric_gas_turb, Re_crit_low, Re_crit_high, Re_gas_part);

// Pure Gas and Pure Liquid Flow Pressure Drop Gradient
// ----------------------------------------------------
   dpdz_liq = -pipeFric_liq/d_hyd / (2*rho_bubble) * m_flow_density_liq * abs(m_flow_density_liq);
   dpdz_gas = -pipeFric_gas/d_hyd / (2*rho_dew)    * m_flow_density_gas * abs(m_flow_density_gas);

// Lockhart-Martinelli Parameter
// -----------------------------
   X2 = if noEvent(X_gas > 0.5) then pipeFric_liq/pipeFric_gas * rho_dew/rho_bubble * (X_liq/X_gas)^2
                                else pipeFric_gas/pipeFric_liq * rho_bubble/rho_dew * (X_gas/X_liq)^2;

//   X2 = pipeFric_liq/pipeFric_gas * rho_dew/rho_bubble * (X_liq/X_gas)^2;

// Weighting Constant of Lockhart-Martinelli Approximation
// -------------------------------------------------------
   C = (5*X_liq_lam + 10*X_liq_turb)*X_liq + (12*X_gas_lam + 20*X_gas_turb)*X_gas;

// Lockhart-Martinelli Approximation
// ---------------------------------
   Phi_2 = 1 + C*sqrt(X2) + X2;

// Friction Pressure Drop Gradient
// -------------------------------
   if noEvent(X_gas > 0.0 and X_gas < 1.0) then

   // 2 Phase Flow Friction Pressure Drop Gradient
   // --------------------------------------------
      dpdz_fric = if noEvent(X_gas > 0.5) then Phi_2 * dpdz_gas
                                          else Phi_2 * dpdz_liq;

   else

   // 1 Phase Flow Friction Pressure Drop Gradient
   // --------------------------------------------
      dpdz_fric = if noEvent(X_gas > 0.5) then -pipeFric_gas/d_hyd / (2*rho) * m_flow_density * abs(m_flow_density)
                                          else -pipeFric_liq/d_hyd / (2*rho) * m_flow_density * abs(m_flow_density);
   end if;

// Friction Pressure Drop
// ----------------------
   dp = dpdz_fric * L_flow;

// Specific Outflow Enthalpy
// -------------------------
   h_out = h;

// No Change of Mass Flow
// ----------------------
   dm_flow = 0;

// No Change of Composition
// ------------------------
   Xi_out = Xi_in;

  annotation (
    Icon(coordinateSystem(preserveAspectRatio=false), graphics={Text(
          extent={{-70,-66},{70,-100}},
          textColor={0,0,0},
          textString="2-Phase Medium")}),
    Diagram(coordinateSystem(preserveAspectRatio=false)),
    Documentation(info="<html>
    For a single conduction element this component models
    <ul>
    <li>the temperature drop or rise of the fluid along the pipe channel</li>
    <li>the heat flow between fluid and pipe wall</li>
    <li>the pressure drop along the pipe channel.</li>
    </ul>
    Condensation and Evaporation are supported.
    
    <p>
    <strong>Media</strong>
    <p>
    The fluid medium is a single substance, which can (partially) change its phase along the pipe flow channel.<br>
    The gas partition can decrease while the liquid partition increases (condensation) or the gas partition can increase while the liquid partition decreases (evaporation).<br>
    Possible fluids are refrigerants like R134a or methanol.<br>
    Furthermore all fluids are allowed which do not change their phase: liquid water, glycol, thermal oil, non-condensing gas.

    <p>
    <strong>Modeling Approach</strong>
    <p>
    <em>Single Phase Flow</em>
    <p>
    The heat exchange rate is calculated according to the Nusselt Number.<br>
    The Nusselt number for flows in pipes is calculated
    <ul>
    <li>in the case of laminar flow  according to Gnielinksi as documented with formula (3.250)</li>
    <li>in the case of turbulent flow according to Gnielinksi as documented with formula (3.260)</li>
    </ul>
    as documented with formula (3.211) in the text book from H.D. Baehr and K. Stephan: \"W&auml;rme- und Stoff&uuml;bertragung\", 7th edition.</li>
    <p>
    The heat exchange rate to the pipe wall is calulated at the beginning of the pipe. Along the length of the pipe an exponential convergence to the wall temperature is assumed.
    
    <p>
    <em>Condensation in Two Phase Flow</em>
    <p>
    Condensation can take place when the fluid is partially in gas phase and the wall temperature is below the dew temperature of the fluid.
    <p>
    The condensation process is modeled according to Nusselt's water film theory. 
    The thickness of the water film around the pipe is determined by the balance of gravity and sheer force.
    The condensation mass flow is determined by the heat exchange through the water film, which depends on the temperature difference between gas and liquid phase.<br>
    
    All the effects are described in chapter 4.1.2 \"Die Nu&szlig;eltsche Wasserhauttheorie\" in the named text book.
    
    <p>
    <em>Evaporation in Two Phase Flow</em>
    <p>
    <p>
    Evaporation can take place when the fluid is partially in liquid phase and the wall temperature is above the boiling temperature of the fluid.
    <p>
    The condensation process is modeled according to Chen's formula (4.160) as described in chapter 4.2.9 \"Zweiphasige Strömungen\" in the named text book.<br>
    There is no physical theory. The model is derived from empirical data. The heat exchange rates depend on the shape of flow and be visualized in flow shape charts.
        
    <p>
    <img src=\"modelica://ThermoFluidStreamPlus/HeatExchangers/HeatExchangersCross/Resources/Evaporation_Process.png\">
    <p>

</html>"),
    __Dymola_Commands);
end ConductionElement_Pipe_2Phase;
