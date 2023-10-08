within TFS_Tutorial.Models;
model SimpleCoolingCycle

  extends Modelica.Icons.Example;

  //Parameters
  parameter Modelica.Units.SI.Temperature heatingThresholdLow=288.15   "Lower threshold for active battery heating" annotation(Dialog(group = "Battery"));
  parameter Modelica.Units.SI.Temperature heatingThresholdHigh=298.15       "Upper threshold for active battery heating" annotation(Dialog(group = "Battery"));
  parameter Modelica.Units.SI.Temperature coolingThreshold=303.15         "Threshold for active cooling" annotation(Dialog(group = "Battery"));
  parameter Modelica.Units.SI.Temperature coolingThresholdVCS=313.15         "Threshold for active VCS cooling" annotation(Dialog(group = "Battery"));
  parameter Modelica.Units.SI.Temperature batteryTargetTemperature=298.15 "Desired battery temperature" annotation(Dialog(group = "Battery"));
  parameter Modelica.Units.SI.Temperature T_amb=288.15 "Ambient temperature" annotation(Dialog(group = "Ambient conditions"));

//Media definitions
  package Coolant =
      ThermofluidStream.Media.myMedia.Incompressible.Examples.Glycol47;
  package Air = ThermofluidStream.Media.myMedia.Air.MoistAir;
  package Refrigerant = ThermofluidStream.Media.XRGMedia.R134a_ph;

  ThermofluidStream.Processes.ConductionElement conductionElement(redeclare package Medium = Coolant, A=10) annotation (Placement(transformation(
        extent={{10,10},{-10,-10}},
        rotation=270,
        origin={-160,0})));
  ThermofluidStream.HeatExchangers.CounterFlowNTU counterFlowNTU(
    redeclare package MediumA = Air,
    redeclare package MediumB = Coolant,
    A=60) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=90,
        origin={60,4})));
  ThermofluidStream.Processes.Pump pump(
    redeclare package Medium = Coolant,
    omega_from_input=true,
    redeclare function dp_tau_pump = ThermofluidStream.Processes.Internal.TurboComponent.dp_tau_centrifugal) annotation (Placement(transformation(extent={{-72,-50},{-92,-30}})));
  ThermofluidStream.Processes.FlowResistance flowResistance(
    redeclare package Medium = Coolant,
    initM_flow=ThermofluidStream.Utilities.Types.InitializationMethods.state,
    l=1,
    shape=ThermofluidStream.Processes.Internal.ShapeOfResistance.circular,
    r(displayUnit="mm") = 0.005,
    redeclare function pLoss = ThermofluidStream.Processes.Internal.FlowResistance.linearQuadraticPressureLoss (k=1000)) annotation (Placement(transformation(extent={{8,58},{28,78}})));
  ThermofluidStream.Boundaries.VolumeFlex volumeFlex(
    redeclare package Medium = Coolant,
    p_start=200000,
    T_start=T_amb) annotation (Placement(transformation(extent={{-106,58},{-86,78}})));
  ThermofluidStream.Boundaries.Source source(redeclare package Medium = Air, temperatureFromInput=true) annotation (Placement(transformation(extent={{136,-50},{116,-30}})));
  ThermofluidStream.Boundaries.Sink sink(redeclare package Medium = Air) annotation (Placement(transformation(extent={{140,58},{160,78}})));
  Modelica.Blocks.Sources.RealExpression realExpression1(y=T_amb) annotation (Placement(transformation(extent={{162,-50},{142,-30}})));
  ThermofluidStream.Processes.Fan fan(
    redeclare package Medium = Air,
    initM_flow=ThermofluidStream.Utilities.Types.InitializationMethods.state,
    omega_from_input=true,
    redeclare function dp_tau_fan = ThermofluidStream.Processes.Internal.TurboComponent.dp_tau_const_isentrop (omega_ref=100, eta=0.8)) annotation (Placement(transformation(extent={{74,58},{94,78}})));
  ThermofluidStream.Processes.FlowResistance flowResistance1(
    redeclare package Medium = Air,
    l=1,
    shape=ThermofluidStream.Processes.Internal.ShapeOfResistance.rectangle,
    a=0.6,
    b=0.2,
    redeclare function pLoss = ThermofluidStream.Processes.Internal.FlowResistance.laminarTurbulentPressureLoss) annotation (Placement(transformation(extent={{110,58},{130,78}})));
  Modelica.Blocks.Sources.Ramp           ramp(
    height=1500,
    duration=1,
    offset=0)
    annotation (Placement(transformation(
        extent={{10,-10},{-10,10}},
        rotation=0,
        origin={124,34})));
  Controllers.LimPID PID(
    controllerType=Modelica.Blocks.Types.SimpleController.PI,
    k=-25,
    Ti=50,
    yMax=1000,
    yMin=5,
    y_start=5) annotation (Placement(transformation(extent={{-46,-82},{-66,-62}})));
  Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor temperatureSensor annotation (Placement(transformation(extent={{-182,-74},{-162,-54}})));
  Modelica.Blocks.Sources.RealExpression realExpression2(y=batteryTargetTemperature) annotation (Placement(transformation(extent={{-14,-82},{-34,-62}})));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow annotation (Placement(transformation(extent={{10,-10},{-10,10}},
        rotation=90,
        origin={-198,24})));
  Modelica.Blocks.Sources.Constant const(k=5000) annotation (Placement(transformation(extent={{-244,44},{-224,64}})));
  ThermofluidStream.Sensors.MultiSensor_Tpm multiSensor_Tpm(
    redeclare package Medium = Coolant,
    temperatureUnit="degC",
    pressureUnit="bar") annotation (Placement(transformation(extent={{-120,-40},{-140,-20}})));
  ThermofluidStream.Sensors.MultiSensor_Tpm multiSensor_Tpm1(
    redeclare package Medium = Coolant,
    temperatureUnit="degC",
    pressureUnit="bar") annotation (Placement(transformation(
        extent={{10,10},{-10,-10}},
        rotation=270,
        origin={-170,42})));
  ThermofluidStream.Sensors.MultiSensor_Tpm multiSensor_Tpm2(
    redeclare package Medium = Coolant,
    temperatureUnit="degC",
    pressureUnit="bar") annotation (Placement(transformation(
        extent={{10,-10},{-10,10}},
        rotation=90,
        origin={44,44})));
  ThermofluidStream.Sensors.MultiSensor_Tpm multiSensor_Tpm3(
    redeclare package Medium = Coolant,
    temperatureUnit="degC",
    pressureUnit="bar") annotation (Placement(transformation(
        extent={{10,-10},{-10,10}},
        rotation=90,
        origin={44,-22})));
  Submodels.Battery battery(ambientTemperature=T_amb) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=90,
        origin={-236,0})));
equation
  connect(volumeFlex.outlet, flowResistance.inlet) annotation (Line(
      points={{-86,68},{8,68}},
      color={28,108,200},
      thickness=0.5));
  connect(source.T0_var, realExpression1.y) annotation (Line(points={{128,-40},{141,-40}}, color={0,0,127}));
  connect(sink.inlet, flowResistance1.outlet) annotation (Line(
      points={{140,68},{130,68}},
      color={28,108,200},
      thickness=0.5));
  connect(fan.outlet, flowResistance1.inlet) annotation (Line(
      points={{94,68},{100,68},{100,68},{110,68}},
      color={28,108,200},
      thickness=0.5));
  connect(counterFlowNTU.outletA, fan.inlet) annotation (Line(
      points={{66,14},{66,68},{74,68}},
      color={28,108,200},
      thickness=0.5));
  connect(source.outlet, counterFlowNTU.inletA) annotation (Line(
      points={{116,-40},{66,-40},{66,-6}},
      color={28,108,200},
      thickness=0.5));
  connect(ramp.y, fan.omega_input) annotation (Line(points={{113,34},{84,34},{84,58}}, color={0,0,127}));
  connect(temperatureSensor.T, PID.u_m) annotation (Line(points={{-161,-64},{-146,-64},{-146,-92},{-56,-92},{-56,-84}}, color={0,0,127}));
  connect(PID.u_s, realExpression2.y) annotation (Line(points={{-44,-72},{-35,-72}}, color={0,0,127}));
  connect(PID.y, pump.omega_input) annotation (Line(points={{-67,-72},{-82,-72},{-82,-50}}, color={0,0,127}));
  connect(prescribedHeatFlow.Q_flow, const.y) annotation (Line(points={{-198,34},{-198,54},{-223,54}}, color={0,0,127}));
  connect(pump.outlet, multiSensor_Tpm.inlet) annotation (Line(
      points={{-92,-40},{-120,-40}},
      color={28,108,200},
      thickness=0.5));
  connect(conductionElement.inlet, multiSensor_Tpm.outlet) annotation (Line(
      points={{-160,-10},{-160,-40},{-140,-40}},
      color={28,108,200},
      thickness=0.5));
  connect(conductionElement.outlet, multiSensor_Tpm1.inlet) annotation (Line(
      points={{-160,10},{-160,32}},
      color={28,108,200},
      thickness=0.5));
  connect(multiSensor_Tpm1.outlet, volumeFlex.inlet) annotation (Line(
      points={{-160,52},{-160,68},{-106,68}},
      color={28,108,200},
      thickness=0.5));
  connect(flowResistance.outlet, multiSensor_Tpm2.inlet) annotation (Line(
      points={{28,68},{54,68},{54,54}},
      color={28,108,200},
      thickness=0.5));
  connect(multiSensor_Tpm2.outlet, counterFlowNTU.inletB) annotation (Line(
      points={{54,34},{54,14}},
      color={28,108,200},
      thickness=0.5));
  connect(counterFlowNTU.outletB, multiSensor_Tpm3.inlet) annotation (Line(
      points={{54,-6},{54,-12}},
      color={28,108,200},
      thickness=0.5));
  connect(multiSensor_Tpm3.outlet, pump.inlet) annotation (Line(
      points={{54,-32},{54,-40},{-72,-40}},
      color={28,108,200},
      thickness=0.5));
  connect(battery.port_a, conductionElement.heatPort) annotation (Line(points={{-226,0},{-169.8,0}}, color={191,0,0}));
  connect(prescribedHeatFlow.port, battery.port_a) annotation (Line(points={{-198,14},{-200,14},{-200,0},{-226,0}}, color={191,0,0}));
  connect(temperatureSensor.port, battery.port_a) annotation (Line(points={{-182,-64},{-176,-64},{-176,-62},{-200,-62},{-200,0},{-226,0}}, color={191,0,0}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-260,-160},{260,160}})),
    experiment(StopTime=3600, __Dymola_Algorithm="Dassl"));
end SimpleCoolingCycle;
