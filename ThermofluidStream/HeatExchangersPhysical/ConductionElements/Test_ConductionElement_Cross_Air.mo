within ThermofluidStream.HeatExchangersPhysical.ConductionElements;
model Test_ConductionElement_Cross_Air "Test of heat exchanger Cross Flow Conduction Element for Non-Condensing Gas"
    extends Modelica.Icons.Example;

  replaceable package Medium_CrossChannel =
      ThermofluidStream.Media.myMedia.Air.MoistAir                        constrainedby ThermofluidStream.Media.myMedia.Interfaces.PartialMedium
                                                        annotation(choicesAllMatching = true, Dialog(group = "Medium definitions"));

  ConductionElements.ConductionElement_Cross_1Phase conductionElement_Cross(redeclare package Medium = Medium_CrossChannel) annotation (Placement(transformation(extent={{-20,-40},{20,0}})));
  ThermofluidStream.Boundaries.Source source_CrossChannel(
    redeclare package Medium = Medium_CrossChannel,
    temperatureFromInput=true,
    xiFromInput=true,
    T0_par=293.15,
    p0_par=100010) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={-132,-20})));
  ThermofluidStream.Boundaries.Sink sink_CrossChannel(redeclare package Medium = Medium_CrossChannel, p0_par=100000) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={140,-20})));
  inner ThermofluidStream.DropOfCommons dropOfCommons(m_flow_reg=1.0E-20,
                                                      assertionLevel=AssertionLevel.warning)
                                                                                annotation (Placement(transformation(extent={{-180,-100},{-160,-80}})));
  ThermofluidStream.FlowControl.MCV mCV_cross(
    redeclare package Medium = Medium_CrossChannel,
    initM_flow=ThermofluidStream.Utilities.Types.InitializationMethods.state,
    m_flow_0=0,
    mode=ThermofluidStream.FlowControl.Internal.Types.MassflowControlValveMode.volume_flow,
    setpointFromInput=true,
    volumeFlow_set_par=0.08333) annotation (Placement(transformation(extent={{-116,-30},{-96,-10}})));
  Modelica.Blocks.Sources.Constant CrossFlowRate(k(
      quantity="VolumeFlowRate",
      unit="m3/s",
      displayUnit="m3/h") = 0.0083333333333333, y(
      quantity="VolumeFlowRate",
      unit="m3/s",
      displayUnit="m3/h")) annotation (Placement(transformation(extent={{-140,20},{-120,40}})));
  Modelica.Blocks.Sources.Constant CrossInflowTemperature(k(
      quantity="Temperature",
      unit="K",
      displayUnit="degC") = 296.35, y(
      quantity="Temperature",
      unit="K",
      displayUnit="degC")) annotation (Placement(transformation(extent={{-180,-20},{-160,0}})));
  Modelica.Blocks.Sources.Constant CrossInflowPressure(k(
      quantity="AbsolutePressure",
      unit="Pa",
      displayUnit="bar") = 200000, y(
      quantity="AbsolutePressure",
      unit="Pa",
      displayUnit="bar")) annotation (Placement(transformation(extent={{-140,60},{-120,80}})));
  Modelica.Blocks.Sources.Constant CrossInflowHumidity(k(
      quantity="MassFraction",
      unit="kg/kg",
      displayUnit="g/kg") = 0.015,y(
      quantity="MassFraction",
      unit="kg/kg",
      displayUnit="g/kg")) annotation (Placement(transformation(extent={{-180,-50},{-160,-30}})));
  ThermofluidStream.Sensors.SingleSensorSelect Enthalpy_AirOut_Sensor(
    redeclare package Medium = Medium_CrossChannel,
    digits=0,
    quantity=ThermofluidStream.Sensors.Internal.Types.Quantities.h_Jpkg) annotation (Placement(transformation(extent={{76,2},{102,22}})));
  ThermofluidStream.Sensors.SingleSensorSelect Temperature_AirOut_Sensor(redeclare package Medium = Medium_CrossChannel, quantity=ThermofluidStream.Sensors.Internal.Types.Quantities.T_C) annotation (Placement(transformation(extent={{78,-10},{98,10}})));
  ThermofluidStream.Sensors.SingleSensorSelect Pressure_CrossIn_Sensor(
    redeclare package Medium = Medium_CrossChannel,
    digits=1,
    quantity=ThermofluidStream.Sensors.Internal.Types.Quantities.p_Pa) annotation (Placement(transformation(extent={{-82,18},{-46,38}})));
  ThermofluidStream.Sensors.SingleSensorSelect Enthalpy_CrossIn_Sensor(
    redeclare package Medium = Medium_CrossChannel,
    digits=0,
    quantity=ThermofluidStream.Sensors.Internal.Types.Quantities.h_Jpkg) annotation (Placement(transformation(extent={{-82,-10},{-56,10}})));
  Modelica.Thermal.HeatTransfer.Celsius.PrescribedTemperature tinTemperature annotation (Placement(transformation(extent={{40,60},{20,80}})));
  Modelica.Blocks.Sources.Constant T_wall(k(
      quantity="Temperature_degC",
      unit="degC") = 60,  y(
      quantity="Temperature_degC",
      unit="degC")) annotation (Placement(transformation(extent={{80,60},{60,80}})));
  HXutilities.Sensors.SingleSensorHumidity phi_CrossIn_Sensor(redeclare package Medium = Medium_CrossChannel, digits=2) annotation (Placement(transformation(extent={{-80,4},{-58,24}})));
  HXutilities.Sensors.SingleSensorHumidity phi_AirOut_Sensor(redeclare package Medium = Medium_CrossChannel, digits=2) annotation (Placement(transformation(extent={{78,16},{98,36}})));
  ThermofluidStream.Sensors.SingleSensorX Humidity_AirOut_SensorX(redeclare package Medium = Medium_CrossChannel, digits=4) annotation (Placement(transformation(extent={{74,30},{114,50}})));
  ThermofluidStream.Sensors.SingleFlowSensor EnthalpyFlow_AirOut_Sensor(
    redeclare package Medium = Medium_CrossChannel,
    digits=1,
    quantity=ThermofluidStream.Sensors.Internal.Types.MassFlowQuantities.H_flow_Jps) annotation (Placement(transformation(extent={{80,-24},{100,-4}})));
  ThermofluidStream.Sensors.SingleFlowSensor MassFlow_AirOut_Sensor(
    redeclare package Medium = Medium_CrossChannel,
    digits=3,
    quantity=ThermofluidStream.Sensors.Internal.Types.MassFlowQuantities.m_flow_gps) annotation (Placement(transformation(extent={{34,-24},{54,-4}})));
  ThermofluidStream.Sensors.SingleFlowSensor EnthalpyFlow_AirIn_Sensor(
    redeclare package Medium = Medium_CrossChannel,
    digits=1,
    quantity=ThermofluidStream.Sensors.Internal.Types.MassFlowQuantities.H_flow_Jps) annotation (Placement(transformation(extent={{-54,-24},{-34,-4}})));
  ThermofluidStream.Sensors.SingleFlowSensor MassFlow_AirIn_Sensor(
    redeclare package Medium = Medium_CrossChannel,
    digits=3,
    quantity=ThermofluidStream.Sensors.Internal.Types.MassFlowQuantities.m_flow_gps) annotation (Placement(transformation(extent={{-82,-24},{-62,-4}})));
  Modelica.Thermal.HeatTransfer.Sensors.HeatFlowSensor heatFlowSensor annotation (Placement(transformation(
        extent={{10,-10},{-10,10}},
        rotation=90,
        origin={0,30})));
  Modelica.Blocks.Interaction.Show.RealValue heatFlowValue(significantDigits=4) annotation (Placement(transformation(extent={{20,20},{40,40}})));
  ThermofluidStream.Sensors.DifferenceSensorSelect pressureDropSensor(
    redeclare package MediumA = Medium_CrossChannel,
    redeclare package MediumB = Medium_CrossChannel,
    digits=3,
    quantity=ThermofluidStream.Sensors.Internal.Types.Quantities.p_Pa) annotation (Placement(transformation(extent={{36,-60},{66,-40}})));
equation
  connect(CrossFlowRate.y, mCV_cross.setpoint_var) annotation (Line(points={{-119,30},{-106,30},{-106,-12}},
                                                                                                          color={0,0,127}));
  connect(source_CrossChannel.outlet, mCV_cross.inlet) annotation (Line(
      points={{-122,-20},{-116,-20}},
      color={28,108,200},
      thickness=0.5));
  connect(source_CrossChannel.T0_var, CrossInflowTemperature.y) annotation (Line(points={{-134,-20},{-148,-20},{-148,-10},{-159,-10}},
                                                                                                                                   color={0,0,127}));
  connect(Pressure_CrossIn_Sensor.inlet, mCV_cross.outlet) annotation (Line(
      points={{-82,28},{-90,28},{-90,-20},{-96,-20}},
      color={28,108,200},
      thickness=0.5));
  connect(Enthalpy_CrossIn_Sensor.inlet, mCV_cross.outlet) annotation (Line(
      points={{-82,0},{-90,0},{-90,-20},{-96,-20}},
      color={28,108,200},
      thickness=0.5));
  connect(CrossInflowHumidity.y, source_CrossChannel.xi_var[1]) annotation (Line(points={{-159,-40},{-148,-40},{-148,-26},{-134,-26}},
                                                                                                                color={0,0,127}));
  connect(tinTemperature.T, T_wall.y) annotation (Line(points={{42,70},{59,70}}, color={0,0,127}));
  connect(mCV_cross.outlet, phi_CrossIn_Sensor.inlet) annotation (Line(
      points={{-96,-20},{-90,-20},{-90,14},{-80,14}},
      color={28,108,200},
      thickness=0.5));
  connect(conductionElement_Cross.outlet, MassFlow_AirOut_Sensor.inlet) annotation (Line(
      points={{20,-20},{34,-20}},
      color={28,108,200},
      thickness=0.5));
  connect(MassFlow_AirOut_Sensor.outlet, EnthalpyFlow_AirOut_Sensor.inlet) annotation (Line(
      points={{54,-20},{80,-20}},
      color={28,108,200},
      thickness=0.5));
  connect(EnthalpyFlow_AirOut_Sensor.outlet, sink_CrossChannel.inlet) annotation (Line(
      points={{100,-20},{130,-20}},
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
  connect(EnthalpyFlow_AirIn_Sensor.outlet, conductionElement_Cross.inlet) annotation (Line(
      points={{-34,-20},{-20,-20}},
      color={28,108,200},
      thickness=0.5));
  connect(MassFlow_AirIn_Sensor.outlet, EnthalpyFlow_AirIn_Sensor.inlet) annotation (Line(
      points={{-62,-20},{-54,-20}},
      color={28,108,200},
      thickness=0.5));
  connect(MassFlow_AirIn_Sensor.inlet, mCV_cross.outlet) annotation (Line(
      points={{-82,-20},{-96,-20}},
      color={28,108,200},
      thickness=0.5));
  connect(heatFlowValue.numberPort, heatFlowSensor.Q_flow) annotation (Line(points={{18.5,30},{11,30}}, color={0,0,127}));
  connect(heatFlowSensor.port_a, tinTemperature.port) annotation (Line(points={{0,40},{0,70},{20,70}}, color={191,0,0}));
  connect(heatFlowSensor.port_b, conductionElement_Cross.heatPort) annotation (Line(points={{0,20},{0,-0.4}}, color={191,0,0}));
  connect(pressureDropSensor.inletA, conductionElement_Cross.outlet) annotation (Line(
      points={{36,-46},{28,-46},{28,-20},{20,-20}},
      color={28,108,200},
      thickness=0.5));
  connect(pressureDropSensor.inletB, EnthalpyFlow_AirIn_Sensor.outlet) annotation (Line(
      points={{36,-54},{-28,-54},{-28,-20},{-34,-20}},
      color={28,108,200},
      thickness=0.5));
  annotation (
    Icon(coordinateSystem(preserveAspectRatio=false, extent={{-180,-100},{140,100}})),
    Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-180,-100},{140,100}}), graphics={
        Text(
          extent={{80,-40},{120,-60}},
          textColor={0,0,0},
          textString="Cross Flow"),
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
          extent={{-36,-10},{-32,-16}},
          textColor={0,0,0},
          textString="W"),
        Text(
          extent={{-142,0},{-118,-14}},
          textColor={0,0,0},
          textString="Moist Air"),
        Text(
          extent={{130,0},{146,-12}},
          textColor={0,0,0},
          textString="Moist Air")}),
    experiment(
      StopTime=100,
      Interval=0.01,
      __Dymola_Algorithm="Esdirk23a"),
    __Dymola_Commands);
end Test_ConductionElement_Cross_Air;
