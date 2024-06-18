within ThermofluidStream.HeatExchangersPhysical.ConductionElements;
model Test_ConductionElement_Pipe_CondensingAir_Sequence "Test of 2 heat exchanger Conduction Elements Pipe for Condensing Gas in a sequence"
    extends Modelica.Icons.Example;

  replaceable package Medium_AirChannel =
      ThermofluidStream.Media.myMedia.Air.MoistAir                        constrainedby ThermofluidStream.Media.myMedia.Interfaces.PartialMedium
                                                        annotation(choicesAllMatching = true, Dialog(group = "Medium definitions"));
  replaceable package Medium_WaterChannel =
      ThermofluidStream.Media.myMedia.Water.StandardWater "Medium model for condensed vapor in cross channel (A)"
    annotation (choicesAllMatching=true, Dialog(group = "Medium definitions"));

  ConductionElements.ConductionElement_Pipe_CondensingAir conductionElement(initM_flow=ThermofluidStream.Utilities.Types.InitializationMethods.state, redeclare package Medium_Water = Medium_WaterChannel) annotation (Placement(transformation(extent={{-20,-40},{20,0}})));
  ThermofluidStream.Boundaries.Source source_AirChannel(
    redeclare package Medium = Medium_AirChannel,
    temperatureFromInput=true,
    xiFromInput=true,
    T0_par=293.15,
    p0_par=100010) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={-132,-20})));
  ThermofluidStream.Boundaries.Sink sink_AirChannel(redeclare package Medium = Medium_AirChannel, p0_par=100000) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={276,-20})));
  inner ThermofluidStream.DropOfCommons dropOfCommons(m_flow_reg=1.0E-20,
                                                      assertionLevel=AssertionLevel.warning)
                                                                                annotation (Placement(transformation(extent={{-180,-100},{-160,-80}})));
  ThermofluidStream.FlowControl.MCV mCV(
    redeclare package Medium = Medium_AirChannel,
    initM_flow=ThermofluidStream.Utilities.Types.InitializationMethods.state,
    m_flow_0=0,
    mode=ThermofluidStream.FlowControl.Internal.Types.MassflowControlValveMode.volume_flow,
    setpointFromInput=true,
    volumeFlow_set_par=0.08333) annotation (Placement(transformation(extent={{-116,-30},{-96,-10}})));
  Modelica.Blocks.Sources.Constant AirFlowRate(k(
      quantity="VolumeFlowRate",
      unit="m3/s",
      displayUnit="m3/h") = 0.0083333333333333, y(
      quantity="VolumeFlowRate",
      unit="m3/s",
      displayUnit="m3/h")) annotation (Placement(transformation(extent={{-140,20},{-120,40}})));
  Modelica.Blocks.Sources.Constant AirInflowTemperature(k(
      quantity="Temperature",
      unit="K",
      displayUnit="degC") = 296.35, y(
      quantity="Temperature",
      unit="K",
      displayUnit="degC")) annotation (Placement(transformation(extent={{-180,-20},{-160,0}})));
  Modelica.Blocks.Sources.Constant AirInflowPressure(k(
      quantity="AbsolutePressure",
      unit="Pa",
      displayUnit="bar") = 200000, y(
      quantity="AbsolutePressure",
      unit="Pa",
      displayUnit="bar")) annotation (Placement(transformation(extent={{-140,60},{-120,80}})));
  Modelica.Blocks.Sources.Constant AirInflowHumidity(k(
      quantity="MassFraction",
      unit="kg/kg",
      displayUnit="g/kg") = 0.015, y(
      quantity="MassFraction",
      unit="kg/kg",
      displayUnit="g/kg")) annotation (Placement(transformation(extent={{-180,-50},{-160,-30}})));
  ThermofluidStream.Sensors.SingleSensorSelect Enthalpy_AirOut_Sensor(
    redeclare package Medium = Medium_AirChannel,
    digits=0,
    quantity=ThermofluidStream.Sensors.Internal.Types.Quantities.h_Jpkg) annotation (Placement(transformation(extent={{76,2},{102,22}})));
  ThermofluidStream.Sensors.SingleSensorSelect Temperature_AirOut_Sensor(redeclare package Medium = Medium_AirChannel, quantity=ThermofluidStream.Sensors.Internal.Types.Quantities.T_C) annotation (Placement(transformation(extent={{78,-10},{98,10}})));
  ThermofluidStream.Sensors.SingleSensorSelect Pressure_AirIn_Sensor(
    redeclare package Medium = Medium_AirChannel,
    digits=1,
    quantity=ThermofluidStream.Sensors.Internal.Types.Quantities.p_Pa) annotation (Placement(transformation(extent={{-82,18},{-46,38}})));
  ThermofluidStream.Sensors.SingleSensorSelect Enthalpy_AirIn_Sensor(
    redeclare package Medium = Medium_AirChannel,
    digits=0,
    quantity=ThermofluidStream.Sensors.Internal.Types.Quantities.h_Jpkg) annotation (Placement(transformation(extent={{-82,-10},{-56,10}})));
  Modelica.Thermal.HeatTransfer.Celsius.PrescribedTemperature tinTemperature annotation (Placement(transformation(extent={{180,60},{160,80}})));
  Modelica.Blocks.Sources.Constant T_wall(k(
      quantity="Temperature_degC",
      unit="degC") = 5.0, y(
      quantity="Temperature_degC",
      unit="degC")) annotation (Placement(transformation(extent={{220,60},{200,80}})));
  ThermofluidStream.Boundaries.TerminalSource terminalSource(redeclare package Medium = Medium_WaterChannel,                                        p_0=100100) annotation (Placement(transformation(extent={{-54,-42},{-34,-22}})));
  HXutilities.Sensors.SingleSensorHumidity phi_AirIn_Sensor(redeclare package Medium = Medium_AirChannel, digits=2) annotation (Placement(transformation(extent={{-80,4},{-58,24}})));
  ThermofluidStream.Boundaries.Sink sink_waterChannel(redeclare package Medium = Medium_WaterChannel, p0_par=100000) annotation (Placement(transformation(extent={{256,-60},{276,-40}})));
  HXutilities.Sensors.SingleSensorHumidity phi_AirOut_Sensor(redeclare package Medium = Medium_AirChannel, digits=2) annotation (Placement(transformation(extent={{78,16},{98,36}})));
  ThermofluidStream.Sensors.SingleSensorX Humidity_AirOut_SensorX(redeclare package Medium = Medium_AirChannel, digits=4) annotation (Placement(transformation(extent={{74,30},{114,50}})));
  ThermofluidStream.Sensors.SingleFlowSensor MassFlow_WaterOut_Sensor(
    redeclare package Medium = Medium_WaterChannel,
    digits=3,
    quantity=ThermofluidStream.Sensors.Internal.Types.MassFlowQuantities.m_flow_gps) annotation (Placement(transformation(extent={{50,-54},{70,-34}})));
  ThermofluidStream.Sensors.SingleFlowSensor EnthalpyFlow_AirOut_Sensor(
    redeclare package Medium = Medium_AirChannel,
    digits=1,
    quantity=ThermofluidStream.Sensors.Internal.Types.MassFlowQuantities.H_flow_Jps) annotation (Placement(transformation(extent={{80,-24},{100,-4}})));
  ThermofluidStream.Sensors.SingleFlowSensor EnthalpyFlow_WaterOut_Sensor(
    redeclare package Medium = Medium_WaterChannel,
    digits=1,
    quantity=ThermofluidStream.Sensors.Internal.Types.MassFlowQuantities.H_flow_Jps) annotation (Placement(transformation(extent={{86,-54},{106,-34}})));
  ThermofluidStream.Sensors.SingleFlowSensor MassFlow_AirOut_Sensor(
    redeclare package Medium = Medium_AirChannel,
    digits=3,
    quantity=ThermofluidStream.Sensors.Internal.Types.MassFlowQuantities.m_flow_gps) annotation (Placement(transformation(extent={{34,-24},{54,-4}})));
  ThermofluidStream.Sensors.SingleFlowSensor EnthalpyFlow_AirIn_Sensor(
    redeclare package Medium = Medium_AirChannel,
    digits=1,
    quantity=ThermofluidStream.Sensors.Internal.Types.MassFlowQuantities.H_flow_Jps) annotation (Placement(transformation(extent={{-54,-24},{-34,-4}})));
  ThermofluidStream.Sensors.SingleFlowSensor MassFlow_AirIn_Sensor(
    redeclare package Medium = Medium_AirChannel,
    digits=3,
    quantity=ThermofluidStream.Sensors.Internal.Types.MassFlowQuantities.m_flow_gps) annotation (Placement(transformation(extent={{-82,-24},{-62,-4}})));
  Modelica.Thermal.HeatTransfer.Sensors.HeatFlowSensor heatFlowSensor annotation (Placement(transformation(
        extent={{10,-10},{-10,10}},
        rotation=90,
        origin={0,30})));
  Modelica.Blocks.Interaction.Show.RealValue heatFlowValue(significantDigits=4) annotation (Placement(transformation(extent={{20,20},{40,40}})));
  ThermofluidStream.Sensors.SingleSensorSelect Temperature_WaterOut_Sensor(redeclare package Medium = Medium_WaterChannel, quantity=ThermofluidStream.Sensors.Internal.Types.Quantities.T_C) annotation (Placement(transformation(extent={{86,-72},{106,-52}})));
  ThermofluidStream.Sensors.DifferenceSensorSelect pressureDropSensor(
    redeclare package MediumA = Medium_AirChannel,
    redeclare package MediumB = Medium_AirChannel,
    digits=3,
    quantity=ThermofluidStream.Sensors.Internal.Types.Quantities.p_Pa) annotation (Placement(transformation(extent={{-2,-80},{28,-60}})));
  ConductionElement_Pipe_CondensingAir conductionElement1(initM_flow=ThermofluidStream.Utilities.Types.InitializationMethods.state, redeclare package Medium_Water = Medium_WaterChannel) annotation (Placement(transformation(extent={{120,-40},{160,0}})));
  ThermofluidStream.Sensors.SingleSensorSelect Enthalpy_AirOut_Sensor1(
    redeclare package Medium = Medium_AirChannel,
    digits=0,
    quantity=ThermofluidStream.Sensors.Internal.Types.Quantities.h_Jpkg) annotation (Placement(transformation(extent={{216,2},{242,22}})));
  ThermofluidStream.Sensors.SingleSensorSelect Temperature_AirOut_Sensor1(redeclare package Medium = Medium_AirChannel, quantity=ThermofluidStream.Sensors.Internal.Types.Quantities.T_C)
                                                                                                                                                                                         annotation (Placement(transformation(extent={{218,-10},{238,10}})));
  HXutilities.Sensors.SingleSensorHumidity phi_AirOut_Sensor1(redeclare package Medium = Medium_AirChannel, digits=2) annotation (Placement(transformation(extent={{218,16},{238,36}})));
  ThermofluidStream.Sensors.SingleSensorX Humidity_AirOut_SensorX1(redeclare package Medium = Medium_AirChannel, digits=4)
                                                                                                                          annotation (Placement(transformation(extent={{214,30},{254,50}})));
  ThermofluidStream.Sensors.SingleFlowSensor MassFlow_WaterOut_Sensor1(
    redeclare package Medium = Medium_WaterChannel,
    digits=3,
    quantity=ThermofluidStream.Sensors.Internal.Types.MassFlowQuantities.m_flow_gps) annotation (Placement(transformation(extent={{190,-54},{210,-34}})));
  ThermofluidStream.Sensors.SingleFlowSensor EnthalpyFlow_AirOut_Sensor1(
    redeclare package Medium = Medium_AirChannel,
    digits=1,
    quantity=ThermofluidStream.Sensors.Internal.Types.MassFlowQuantities.H_flow_Jps) annotation (Placement(transformation(extent={{220,-24},{240,-4}})));
  ThermofluidStream.Sensors.SingleFlowSensor EnthalpyFlow_WaterOut_Sensor1(
    redeclare package Medium = Medium_WaterChannel,
    digits=1,
    quantity=ThermofluidStream.Sensors.Internal.Types.MassFlowQuantities.H_flow_Jps) annotation (Placement(transformation(extent={{226,-54},{246,-34}})));
  ThermofluidStream.Sensors.SingleFlowSensor MassFlow_AirOut_Sensor1(
    redeclare package Medium = Medium_AirChannel,
    digits=3,
    quantity=ThermofluidStream.Sensors.Internal.Types.MassFlowQuantities.m_flow_gps) annotation (Placement(transformation(extent={{174,-24},{194,-4}})));
  Modelica.Thermal.HeatTransfer.Sensors.HeatFlowSensor heatFlowSensor1
                                                                      annotation (Placement(transformation(
        extent={{10,-10},{-10,10}},
        rotation=90,
        origin={140,30})));
  Modelica.Blocks.Interaction.Show.RealValue heatFlowValue1(significantDigits=4)
                                                                                annotation (Placement(transformation(extent={{160,20},{180,40}})));
  ThermofluidStream.Sensors.SingleSensorSelect Temperature_WaterOut_Sensor1(redeclare package Medium = Medium_WaterChannel, quantity=ThermofluidStream.Sensors.Internal.Types.Quantities.T_C)
                                                                                                                                                                                             annotation (Placement(transformation(extent={{226,-72},{246,-52}})));
  ThermofluidStream.Sensors.DifferenceSensorSelect pressureDropSensor1(
    redeclare package MediumA = Medium_AirChannel,
    redeclare package MediumB = Medium_AirChannel,
    digits=3,
    quantity=ThermofluidStream.Sensors.Internal.Types.Quantities.p_Pa) annotation (Placement(transformation(extent={{138,-80},{168,-60}})));
equation
  connect(AirFlowRate.y, mCV.setpoint_var) annotation (Line(points={{-119,30},{-106,30},{-106,-12}}, color={0,0,127}));
  connect(source_AirChannel.outlet, mCV.inlet) annotation (Line(
      points={{-122,-20},{-116,-20}},
      color={28,108,200},
      thickness=0.5));
  connect(source_AirChannel.T0_var, AirInflowTemperature.y) annotation (Line(points={{-134,-20},{-148,-20},{-148,-10},{-159,-10}}, color={0,0,127}));
  connect(Pressure_AirIn_Sensor.inlet, mCV.outlet) annotation (Line(
      points={{-82,28},{-90,28},{-90,-20},{-96,-20}},
      color={28,108,200},
      thickness=0.5));
  connect(Enthalpy_AirIn_Sensor.inlet, mCV.outlet) annotation (Line(
      points={{-82,0},{-90,0},{-90,-20},{-96,-20}},
      color={28,108,200},
      thickness=0.5));
  connect(AirInflowHumidity.y, source_AirChannel.xi_var[1]) annotation (Line(points={{-159,-40},{-148,-40},{-148,-26},{-134,-26}}, color={0,0,127}));
  connect(tinTemperature.T, T_wall.y) annotation (Line(points={{182,70},{199,70}},
                                                                                 color={0,0,127}));
  connect(terminalSource.outlet, conductionElement.inlet_water) annotation (Line(
      points={{-34,-32},{-20,-32}},
      color={28,108,200},
      thickness=0.5));
  connect(mCV.outlet, phi_AirIn_Sensor.inlet) annotation (Line(
      points={{-96,-20},{-90,-20},{-90,14},{-80,14}},
      color={28,108,200},
      thickness=0.5));
  connect(conductionElement.outlet_water, MassFlow_WaterOut_Sensor.inlet) annotation (Line(
      points={{20,-32},{34,-32},{34,-50},{50,-50}},
      color={28,108,200},
      thickness=0.5));
  connect(MassFlow_WaterOut_Sensor.outlet, EnthalpyFlow_WaterOut_Sensor.inlet) annotation (Line(
      points={{70,-50},{86,-50}},
      color={28,108,200},
      thickness=0.5));
  connect(conductionElement.outlet, MassFlow_AirOut_Sensor.inlet) annotation (Line(
      points={{20,-20},{34,-20}},
      color={28,108,200},
      thickness=0.5));
  connect(MassFlow_AirOut_Sensor.outlet, EnthalpyFlow_AirOut_Sensor.inlet) annotation (Line(
      points={{54,-20},{80,-20}},
      color={28,108,200},
      thickness=0.5));
  connect(Temperature_AirOut_Sensor.inlet, MassFlow_AirOut_Sensor.outlet) annotation (Line(
      points={{78,0},{64,0},{64,-20},{54,-20}},
      color={28,108,200},
      thickness=0.5));
  connect(Enthalpy_AirOut_Sensor.inlet, MassFlow_AirOut_Sensor.outlet) annotation (Line(
      points={{76,12},{64,12},{64,-20},{54,-20}},
      color={28,108,200},
      thickness=0.5));
  connect(phi_AirOut_Sensor.inlet, MassFlow_AirOut_Sensor.outlet) annotation (Line(
      points={{78,26},{64,26},{64,-20},{54,-20}},
      color={28,108,200},
      thickness=0.5));
  connect(Humidity_AirOut_SensorX.inlet, MassFlow_AirOut_Sensor.outlet) annotation (Line(
      points={{74,40},{64,40},{64,-20},{54,-20}},
      color={28,108,200},
      thickness=0.5));
  connect(EnthalpyFlow_AirIn_Sensor.outlet, conductionElement.inlet) annotation (Line(
      points={{-34,-20},{-20,-20}},
      color={28,108,200},
      thickness=0.5));
  connect(MassFlow_AirIn_Sensor.outlet, EnthalpyFlow_AirIn_Sensor.inlet) annotation (Line(
      points={{-62,-20},{-54,-20}},
      color={28,108,200},
      thickness=0.5));
  connect(MassFlow_AirIn_Sensor.inlet, mCV.outlet) annotation (Line(
      points={{-82,-20},{-96,-20}},
      color={28,108,200},
      thickness=0.5));
  connect(heatFlowValue.numberPort, heatFlowSensor.Q_flow) annotation (Line(points={{18.5,30},{11,30}}, color={0,0,127}));
  connect(Temperature_WaterOut_Sensor.inlet, MassFlow_WaterOut_Sensor.outlet) annotation (Line(
      points={{86,-62},{78,-62},{78,-50},{70,-50}},
      color={28,108,200},
      thickness=0.5));
  connect(heatFlowSensor.port_a, tinTemperature.port) annotation (Line(points={{0,40},{0,70},{160,70}},color={191,0,0}));
  connect(heatFlowSensor.port_b, conductionElement.heatPort) annotation (Line(points={{0,20},{0,-0.4}}, color={191,0,0}));
  connect(pressureDropSensor.inletA, conductionElement.outlet) annotation (Line(
      points={{-2,-66},{-10,-66},{-10,-44},{26,-44},{26,-20},{20,-20}},
      color={28,108,200},
      thickness=0.5));
  connect(pressureDropSensor.inletB, EnthalpyFlow_AirIn_Sensor.outlet) annotation (Line(
      points={{-2,-74},{-28,-74},{-28,-20},{-34,-20}},
      color={28,108,200},
      thickness=0.5));
  connect(conductionElement1.outlet_water, MassFlow_WaterOut_Sensor1.inlet) annotation (Line(
      points={{160,-32},{174,-32},{174,-50},{190,-50}},
      color={28,108,200},
      thickness=0.5));
  connect(MassFlow_WaterOut_Sensor1.outlet, EnthalpyFlow_WaterOut_Sensor1.inlet) annotation (Line(
      points={{210,-50},{226,-50}},
      color={28,108,200},
      thickness=0.5));
  connect(conductionElement1.outlet, MassFlow_AirOut_Sensor1.inlet) annotation (Line(
      points={{160,-20},{174,-20}},
      color={28,108,200},
      thickness=0.5));
  connect(MassFlow_AirOut_Sensor1.outlet, EnthalpyFlow_AirOut_Sensor1.inlet) annotation (Line(
      points={{194,-20},{220,-20}},
      color={28,108,200},
      thickness=0.5));
  connect(Temperature_AirOut_Sensor1.inlet, MassFlow_AirOut_Sensor1.outlet) annotation (Line(
      points={{218,0},{204,0},{204,-20},{194,-20}},
      color={28,108,200},
      thickness=0.5));
  connect(Enthalpy_AirOut_Sensor1.inlet, MassFlow_AirOut_Sensor1.outlet) annotation (Line(
      points={{216,12},{204,12},{204,-20},{194,-20}},
      color={28,108,200},
      thickness=0.5));
  connect(phi_AirOut_Sensor1.inlet, MassFlow_AirOut_Sensor1.outlet) annotation (Line(
      points={{218,26},{204,26},{204,-20},{194,-20}},
      color={28,108,200},
      thickness=0.5));
  connect(Humidity_AirOut_SensorX1.inlet, MassFlow_AirOut_Sensor1.outlet) annotation (Line(
      points={{214,40},{204,40},{204,-20},{194,-20}},
      color={28,108,200},
      thickness=0.5));
  connect(heatFlowValue1.numberPort, heatFlowSensor1.Q_flow) annotation (Line(points={{158.5,30},{151,30}}, color={0,0,127}));
  connect(Temperature_WaterOut_Sensor1.inlet, MassFlow_WaterOut_Sensor1.outlet) annotation (Line(
      points={{226,-62},{218,-62},{218,-50},{210,-50}},
      color={28,108,200},
      thickness=0.5));
  connect(heatFlowSensor1.port_a, tinTemperature.port) annotation (Line(points={{140,40},{140,70},{160,70}}, color={191,0,0}));
  connect(heatFlowSensor1.port_b, conductionElement1.heatPort) annotation (Line(points={{140,20},{140,-0.4}}, color={191,0,0}));
  connect(pressureDropSensor1.inletA, conductionElement1.outlet) annotation (Line(
      points={{138,-66},{130,-66},{130,-44},{166,-44},{166,-20},{160,-20}},
      color={28,108,200},
      thickness=0.5));
  connect(EnthalpyFlow_AirOut_Sensor.outlet, conductionElement1.inlet) annotation (Line(
      points={{100,-20},{120,-20}},
      color={28,108,200},
      thickness=0.5));
  connect(EnthalpyFlow_WaterOut_Sensor.outlet, conductionElement1.inlet_water) annotation (Line(
      points={{106,-50},{112,-50},{112,-32},{120,-32}},
      color={28,108,200},
      thickness=0.5));
  connect(EnthalpyFlow_AirOut_Sensor1.outlet, sink_AirChannel.inlet) annotation (Line(
      points={{240,-20},{266,-20}},
      color={28,108,200},
      thickness=0.5));
  connect(EnthalpyFlow_WaterOut_Sensor1.outlet, sink_waterChannel.inlet) annotation (Line(
      points={{246,-50},{256,-50}},
      color={28,108,200},
      thickness=0.5));
  connect(EnthalpyFlow_AirOut_Sensor.outlet, pressureDropSensor1.inletB) annotation (Line(
      points={{100,-20},{110,-20},{110,-74},{138,-74}},
      color={28,108,200},
      thickness=0.5));
  annotation (
    Icon(coordinateSystem(preserveAspectRatio=false, extent={{-180,-100},{280,100}})),
    Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-180,-100},{280,100}}), graphics={
        Text(
          extent={{68,-40},{74,-46}},
          textColor={0,0,0},
          textString="g/s"),
        Text(
          extent={{54,-10},{60,-16}},
          textColor={0,0,0},
          textString="g/s"),
        Text(
          extent={{-64,-10},{-58,-16}},
          textColor={0,0,0},
          textString="g/s"),
        Text(
          extent={{28,24},{34,18}},
          textColor={0,0,0},
          textString="W"),
        Text(
          extent={{104,-40},{108,-46}},
          textColor={0,0,0},
          textString="W"),
        Text(
          extent={{-36,-10},{-32,-16}},
          textColor={0,0,0},
          textString="W"),
        Text(
          extent={{-142,0},{-118,-14}},
          textColor={0,0,0},
          textString="Moist Air"),
        Text(
          extent={{266,0},{282,-12}},
          textColor={0,0,0},
          textString="Moist Air"),
        Text(
          extent={{250,-28},{274,-44}},
          textColor={0,0,0},
          textString="Liquid Water"),
        Text(
          extent={{208,-40},{214,-46}},
          textColor={0,0,0},
          textString="g/s"),
        Text(
          extent={{194,-10},{200,-16}},
          textColor={0,0,0},
          textString="g/s"),
        Text(
          extent={{168,24},{174,18}},
          textColor={0,0,0},
          textString="W"),
        Text(
          extent={{244,-40},{248,-46}},
          textColor={0,0,0},
          textString="W")}),
    experiment(
      StopTime=100,
      Interval=0.01,
      __Dymola_Algorithm="Esdirk23a"),
    __Dymola_Commands);
end Test_ConductionElement_Pipe_CondensingAir_Sequence;
