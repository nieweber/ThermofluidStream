within TFS_Tutorial.Models;
model CoolingCycleWithHeatingAndVCS

  extends Modelica.Icons.Example;

  //Parameters
  parameter Modelica.Units.SI.Temperature heatingThresholdLow=283.15
                                                                "Lower threshold for active battery heating" annotation(Dialog(group = "Battery"));
  parameter Modelica.Units.SI.Temperature heatingThresholdHigh=293.15
                                                                     "Upper threshold for active battery heating" annotation(Dialog(group = "Battery"));
  parameter Modelica.Units.SI.Temperature coolingThreshold=303.15         "Threshold for active cooling" annotation(Dialog(group = "Battery"));
  parameter Modelica.Units.SI.Temperature coolingThresholdVCS=313.15         "Threshold for active VCS cooling" annotation(Dialog(group = "Battery"));
  parameter Modelica.Units.SI.Temperature batteryTargetTemperature=298.15 "Desired battery temperature" annotation(Dialog(group = "Battery"));
  parameter Modelica.Units.SI.Temperature T_amb=308.15 "Ambient temperature" annotation(Dialog(group = "Ambient conditions"));

//Media definitions
  package Coolant =
      ThermofluidStream.Media.myMedia.Incompressible.Examples.Glycol47;
  package Air = ThermofluidStream.Media.myMedia.Air.MoistAir;
  package Refrigerant = ThermofluidStream.Media.XRGMedia.R134a_ph;

  ThermofluidStream.Processes.ConductionElement conductionElement(redeclare package Medium = Coolant, A=10) annotation (Placement(transformation(
        extent={{10,10},{-10,-10}},
        rotation=270,
        origin={-160,20})));
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
    redeclare function dp_tau_pump = ThermofluidStream.Processes.Internal.TurboComponent.dp_tau_centrifugal) annotation (Placement(transformation(extent={{-66,-106},{-86,-86}})));
  ThermofluidStream.Processes.FlowResistance flowResistance(
    redeclare package Medium = Coolant,
    initM_flow=ThermofluidStream.Utilities.Types.InitializationMethods.state,
    l=1,
    shape=ThermofluidStream.Processes.Internal.ShapeOfResistance.circular,
    r(displayUnit="mm") = 0.005,
    redeclare function pLoss = ThermofluidStream.Processes.Internal.FlowResistance.linearQuadraticPressureLoss (k=1000)) annotation (Placement(transformation(extent={{2,120},{22,140}})));
  ThermofluidStream.Boundaries.VolumeFlex volumeFlex(
    redeclare package Medium = Coolant,
    p_start=200000,
    T_start=T_amb) annotation (Placement(transformation(extent={{-112,120},{-92,140}})));
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
  Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor temperatureSensor annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=270,
        origin={-200,-48})));
  Modelica.Blocks.Sources.RealExpression realExpression2(y=batteryTargetTemperature) annotation (Placement(transformation(extent={{18,-130},{-2,-110}})));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow annotation (Placement(transformation(extent={{10,-10},{-10,10}},
        rotation=90,
        origin={-200,74})));
  Modelica.Blocks.Sources.Constant const(k=8000)
                                                annotation (Placement(transformation(extent={{-246,88},{-226,108}})));
  ThermofluidStream.Sensors.MultiSensor_Tpm multiSensor_Tpm(
    redeclare package Medium = Coolant,
    temperatureUnit="degC",
    pressureUnit="bar") annotation (Placement(transformation(extent={{-100,-96},{-120,-76}})));
  ThermofluidStream.Sensors.MultiSensor_Tpm multiSensor_Tpm1(
    redeclare package Medium = Coolant,
    temperatureUnit="degC",
    pressureUnit="bar") annotation (Placement(transformation(
        extent={{10,10},{-10,-10}},
        rotation=270,
        origin={-170,98})));
  ThermofluidStream.Sensors.MultiSensor_Tpm multiSensor_Tpm2(
    redeclare package Medium = Coolant,
    temperatureUnit="degC",
    pressureUnit="bar") annotation (Placement(transformation(
        extent={{10,-10},{-10,10}},
        rotation=90,
        origin={44,48})));
  ThermofluidStream.Sensors.MultiSensor_Tpm multiSensor_Tpm3(
    redeclare package Medium = Coolant,
    temperatureUnit="degC",
    pressureUnit="bar") annotation (Placement(transformation(
        extent={{10,-10},{-10,10}},
        rotation=90,
        origin={44,-28})));
  Submodels.PTC_Heater pTC_Heater(redeclare package Medium = Coolant, T_start=T_amb) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=270,
        origin={-156,-26})));
  ThermofluidStream.FlowControl.Switch switch(redeclare package Medium = Coolant) annotation (Placement(transformation(
        extent={{10,-10},{-10,10}},
        rotation=90,
        origin={54,92})));
  ThermofluidStream.Topology.JunctionT2 junctionT2_1(redeclare package Medium = Coolant) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=90,
        origin={54,-62})));
  ThermofluidStream.Processes.FlowResistance flowResistance2(
    redeclare package Medium = Coolant,
    initM_flow=ThermofluidStream.Utilities.Types.InitializationMethods.state,
    l=1,
    shape=ThermofluidStream.Processes.Internal.ShapeOfResistance.circular,
    r(displayUnit="mm") = 0.005,
    redeclare function pLoss = ThermofluidStream.Processes.Internal.FlowResistance.linearQuadraticPressureLoss (k=1000)) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=270,
        origin={8,20})));
  Modelica.Blocks.Logical.Hysteresis hysteresisHeating(uLow=heatingThresholdLow, uHigh=heatingThresholdHigh) annotation (Placement(transformation(extent={{-110,-36},{-130,-16}})));
  Modelica.Blocks.Logical.Hysteresis hysteresisCooling(uLow=heatingThresholdLow, uHigh=coolingThreshold) annotation (Placement(transformation(extent={{154,96},{134,116}})));
  Modelica.Blocks.Math.BooleanToReal booleanToReal(realTrue=0, realFalse=1) annotation (Placement(transformation(extent={{116,96},{96,116}})));
  Modelica.Blocks.Sources.RealExpression realExpression3(y=temperatureSensor.T)
                                                                  annotation (Placement(transformation(extent={{192,96},{172,116}})));
  Controllers.PumpController pumpController(Ti_PI=10) annotation (Placement(transformation(extent={{-48,-130},{-28,-110}})));
  Modelica.Blocks.Logical.Or or1 annotation (Placement(transformation(extent={{-106,-150},{-86,-130}})));
  Modelica.Blocks.Sources.BooleanExpression booleanExpression(y=not hysteresisHeating.y) annotation (Placement(transformation(extent={{-154,-140},{-134,-120}})));
  Modelica.Blocks.Sources.BooleanExpression booleanExpression1(y=hysteresisCooling.y) annotation (Placement(transformation(extent={{-156,-158},{-136,-138}})));
  Submodels.VCS vCS annotation (Placement(transformation(
        extent={{66,-64},{-66,64}},
        rotation=270,
        origin={276,20})));
  ThermofluidStream.FlowControl.Switch switch1(redeclare package Medium = Coolant)
                                                                                  annotation (Placement(transformation(
        extent={{10,-10},{-10,10}},
        rotation=180,
        origin={54,130})));
  ThermofluidStream.Topology.JunctionT2 junctionT2_2(redeclare package Medium = Coolant) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={54,-96})));
  ThermofluidStream.Processes.FlowResistance flowResistance3(
    redeclare package Medium = Coolant,
    initM_flow=ThermofluidStream.Utilities.Types.InitializationMethods.state,
    l=1,
    shape=ThermofluidStream.Processes.Internal.ShapeOfResistance.circular,
    r(displayUnit="mm") = 0.005,
    redeclare function pLoss = ThermofluidStream.Processes.Internal.FlowResistance.linearQuadraticPressureLoss (k=1000)) annotation (Placement(transformation(extent={{138,-106},{118,-86}})));
  ThermofluidStream.Boundaries.Sink sink1(redeclare package Medium = Air)
                                                                         annotation (Placement(transformation(extent={{400,100},{420,120}})));
  ThermofluidStream.Processes.Fan fan1(
    redeclare package Medium = Air,
    initM_flow=ThermofluidStream.Utilities.Types.InitializationMethods.state,
    omega_from_input=true,
    redeclare function dp_tau_fan = ThermofluidStream.Processes.Internal.TurboComponent.dp_tau_const_isentrop (omega_ref=100, eta=0.8)) annotation (Placement(transformation(extent={{334,120},{354,100}})));
  ThermofluidStream.Processes.FlowResistance flowResistance4(
    redeclare package Medium = Air,
    l=1,
    shape=ThermofluidStream.Processes.Internal.ShapeOfResistance.rectangle,
    a=0.6,
    b=0.2,
    redeclare function pLoss = ThermofluidStream.Processes.Internal.FlowResistance.laminarTurbulentPressureLoss) annotation (Placement(transformation(extent={{370,100},{390,120}})));
  Modelica.Blocks.Sources.Ramp           ramp2(
    height=1500,
    duration=1,
    offset=0)
    annotation (Placement(transformation(
        extent={{10,-10},{-10,10}},
        rotation=0,
        origin={378,134})));
  ThermofluidStream.Boundaries.Source source1(redeclare package Medium = Air, temperatureFromInput=true)
                                                                                                        annotation (Placement(transformation(extent={{382,-80},{362,-60}})));
  Modelica.Blocks.Sources.RealExpression realExpression4(y=T_amb) annotation (Placement(transformation(extent={{408,-80},{388,-60}})));
  ThermofluidStream.Sensors.MultiSensor_Tpm multiSensor_Tpm4(
    redeclare package Medium = Coolant,
    temperatureUnit="degC",
    pressureUnit="bar",
    outputTemperature=true)
                        annotation (Placement(transformation(extent={{10,-10},{-10,10}},
        rotation=90,
        origin={218,-74})));
  Modelica.Blocks.Sources.RealExpression realExpression5(y=10)    annotation (Placement(transformation(extent={{316,98},{296,118}})));
  Modelica.Blocks.Sources.RealExpression realExpression6(y=temperatureSensor.T)
                                                                  annotation (Placement(transformation(extent={{184,138},{164,158}})));
  Modelica.Blocks.Logical.Hysteresis hysteresisCoolingVCS(uLow=heatingThresholdLow, uHigh=coolingThresholdVCS) annotation (Placement(transformation(extent={{146,138},{126,158}})));
  Modelica.Blocks.Math.BooleanToReal booleanToReal1(realTrue=0, realFalse=1) annotation (Placement(transformation(extent={{110,138},{90,158}})));
  Submodels.Battery battery(ambientTemperature=T_amb) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=90,
        origin={-242,20})));
equation
  connect(volumeFlex.outlet, flowResistance.inlet) annotation (Line(
      points={{-92,130},{2,130}},
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
  connect(prescribedHeatFlow.Q_flow,const. y) annotation (Line(points={{-200,84},{-200,98},{-225,98}}, color={0,0,127}));
  connect(pump.outlet, multiSensor_Tpm.inlet) annotation (Line(
      points={{-86,-96},{-100,-96}},
      color={28,108,200},
      thickness=0.5));
  connect(conductionElement.outlet, multiSensor_Tpm1.inlet) annotation (Line(
      points={{-160,30},{-160,88}},
      color={28,108,200},
      thickness=0.5));
  connect(multiSensor_Tpm1.outlet, volumeFlex.inlet) annotation (Line(
      points={{-160,108},{-160,130},{-112,130}},
      color={28,108,200},
      thickness=0.5));
  connect(multiSensor_Tpm2.outlet, counterFlowNTU.inletB) annotation (Line(
      points={{54,38},{54,14}},
      color={28,108,200},
      thickness=0.5));
  connect(counterFlowNTU.outletB, multiSensor_Tpm3.inlet) annotation (Line(
      points={{54,-6},{54,-18}},
      color={28,108,200},
      thickness=0.5));
  connect(conductionElement.inlet, pTC_Heater.outlet) annotation (Line(
      points={{-160,10},{-160,-15.8}},
      color={28,108,200},
      thickness=0.5));
  connect(pTC_Heater.inlet, multiSensor_Tpm.outlet) annotation (Line(
      points={{-160,-35.8},{-160,-96},{-120,-96}},
      color={28,108,200},
      thickness=0.5));
  connect(switch.outletA, multiSensor_Tpm2.inlet) annotation (Line(
      points={{54,82},{54,58}},
      color={28,108,200},
      thickness=0.5));
  connect(multiSensor_Tpm3.outlet, junctionT2_1.inletB) annotation (Line(
      points={{54,-38},{54,-52}},
      color={28,108,200},
      thickness=0.5));
  connect(switch.outletB, flowResistance2.inlet) annotation (Line(
      points={{44,92},{8,92},{8,30}},
      color={28,108,200},
      thickness=0.5));
  connect(flowResistance2.outlet, junctionT2_1.inletA) annotation (Line(
      points={{8,10},{8,-62},{44,-62}},
      color={28,108,200},
      thickness=0.5));
  connect(hysteresisHeating.y, pTC_Heater.u) annotation (Line(points={{-131,-26},{-145.4,-26}}, color={255,0,255}));
  connect(temperatureSensor.T, hysteresisHeating.u) annotation (Line(points={{-200,-59},{-200,-72},{-96,-72},{-96,-26},{-108,-26}}, color={0,0,127}));
  connect(hysteresisCooling.y, booleanToReal.u) annotation (Line(points={{133,106},{118,106}}, color={255,0,255}));
  connect(booleanToReal.y, switch.u) annotation (Line(points={{95,106},{80,106},{80,92},{62,92}}, color={0,0,127}));
  connect(hysteresisCooling.u, realExpression3.y) annotation (Line(points={{156,106},{171,106}}, color={0,0,127}));
  connect(pumpController.y, pump.omega_input) annotation (Line(points={{-48.7143,-120},{-76,-120},{-76,-106}}, color={0,0,127}));
  connect(pumpController.u_s, realExpression2.y) annotation (Line(points={{-28.2857,-120},{-3,-120}}, color={0,0,127}));
  connect(temperatureSensor.T, pumpController.u_m) annotation (Line(points={{-200,-59},{-200,-154},{-38,-154},{-38,-129.714}}, color={0,0,127}));
  connect(or1.y, pumpController.u_switch) annotation (Line(points={{-85,-140},{-45.1429,-140},{-45.1429,-129.714}}, color={255,0,255}));
  connect(booleanExpression1.y, or1.u2) annotation (Line(points={{-135,-148},{-108,-148}}, color={255,0,255}));
  connect(booleanExpression.y, or1.u1) annotation (Line(points={{-133,-130},{-128,-130},{-128,-140},{-108,-140}}, color={255,0,255}));
  connect(flowResistance.outlet, switch1.inlet) annotation (Line(
      points={{22,130},{44,130}},
      color={28,108,200},
      thickness=0.5));
  connect(switch1.outletB, switch.inlet) annotation (Line(
      points={{54,120},{54,102}},
      color={28,108,200},
      thickness=0.5));
  connect(switch1.outletA, vCS.inletAirEvaporator) annotation (Line(
      points={{64,130},{228,130},{228,85.175}},
      color={28,108,200},
      thickness=0.5));
  connect(junctionT2_2.outlet, pump.inlet) annotation (Line(
      points={{44,-96},{-66,-96}},
      color={28,108,200},
      thickness=0.5));
  connect(junctionT2_1.outlet, junctionT2_2.inletA) annotation (Line(
      points={{54,-72},{54,-86},{54,-86}},
      color={28,108,200},
      thickness=0.5));
  connect(flowResistance3.outlet, junctionT2_2.inletB) annotation (Line(
      points={{118,-96},{64,-96}},
      color={28,108,200},
      thickness=0.5));
  connect(sink1.inlet, flowResistance4.outlet) annotation (Line(
      points={{400,110},{390,110}},
      color={28,108,200},
      thickness=0.5));
  connect(fan1.outlet, flowResistance4.inlet) annotation (Line(
      points={{354,110},{370,110}},
      color={28,108,200},
      thickness=0.5));
  connect(ramp2.y, fan1.omega_input) annotation (Line(points={{367,134},{344,134},{344,120}}, color={0,0,127}));
  connect(source1.T0_var, realExpression4.y) annotation (Line(points={{374,-70},{387,-70}}, color={0,0,127}));
  connect(source1.outlet, vCS.inletAirCondenser) annotation (Line(
      points={{362,-70},{324,-70},{324,-45.175}},
      color={28,108,200},
      thickness=0.5));
  connect(vCS.outletAirCondenser, fan1.inlet) annotation (Line(
      points={{324,86},{324,110},{334,110}},
      color={28,108,200},
      thickness=0.5));
  connect(vCS.outletAirEvaporator, multiSensor_Tpm4.inlet) annotation (Line(
      points={{228,-46.825},{228,-64}},
      color={28,108,200},
      thickness=0.5));
  connect(multiSensor_Tpm4.outlet, flowResistance3.inlet) annotation (Line(
      points={{228,-84},{228,-96},{138,-96}},
      color={28,108,200},
      thickness=0.5));
  connect(multiSensor_Tpm4.T_out, vCS.u_m) annotation (Line(points={{212,-84},{212,-88},{196,-88},{196,104},{263.6,104},{263.6,83.9375}}, color={0,0,127}));
  connect(realExpression5.y, vCS.u_set) annotation (Line(points={{295,108},{287.6,108},{287.6,83.9375}}, color={0,0,127}));
  connect(booleanToReal1.u, hysteresisCoolingVCS.y) annotation (Line(points={{112,148},{125,148}}, color={255,0,255}));
  connect(hysteresisCoolingVCS.u, realExpression6.y) annotation (Line(points={{148,148},{163,148}}, color={0,0,127}));
  connect(battery.port_a, conductionElement.heatPort) annotation (Line(points={{-232,20},{-169.8,20}}, color={191,0,0}));
  connect(prescribedHeatFlow.port, battery.port_a) annotation (Line(points={{-200,64},{-200,20},{-232,20}}, color={191,0,0}));
  connect(temperatureSensor.port, battery.port_a) annotation (Line(points={{-200,-38},{-200,20},{-232,20}}, color={191,0,0}));
  connect(booleanToReal1.y, switch1.u) annotation (Line(points={{89,148},{54,148},{54,138}}, color={0,0,127}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-260,-160},{420,160}})));
end CoolingCycleWithHeatingAndVCS;
