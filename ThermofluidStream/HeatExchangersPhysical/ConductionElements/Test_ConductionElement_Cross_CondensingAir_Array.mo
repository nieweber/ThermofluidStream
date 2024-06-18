within ThermofluidStream.HeatExchangersPhysical.ConductionElements;
model Test_ConductionElement_Cross_CondensingAir_Array "Test of heat exchanger cell pair: Condensing Air to Refrigerant R134a"
    extends Modelica.Icons.Example;

  replaceable package Medium_CrossChannel =
      ThermofluidStream.Media.myMedia.Air.MoistAir                        constrainedby ThermofluidStream.Media.myMedia.Interfaces.PartialMedium
                                                        annotation(choicesAllMatching = true);

  replaceable package Medium_Water =
      ThermofluidStream.Media.myMedia.Water.StandardWater "Medium model for condensed vapor in pipe channel"
    annotation (choicesAllMatching=true, Dialog(group = "Medium definitions"));

  parameter Integer nCells = 10;

  ConductionElements.ConductionElement_Cross_CondensingAir conductionElement_Cross[nCells](redeclare package Medium_Water = Medium_Water)
                                                                                           annotation (Placement(transformation(extent={{20,-40},{60,0}})));
  ThermofluidStream.Boundaries.Source source_CrossChannel(
    redeclare package Medium = Medium_CrossChannel,
    temperatureFromInput=true,
    xiFromInput=true,
    T0_par=293.15,
    p0_par=100010) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={-92,-20})));
  ThermofluidStream.Boundaries.Sink sink_CrossChannel(redeclare package Medium = Medium_CrossChannel, p0_par=100000) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={200,-20})));
  inner ThermofluidStream.DropOfCommons dropOfCommons(m_flow_reg=1.0E-20,
                                                      assertionLevel=AssertionLevel.warning)
                                                                                annotation (Placement(transformation(extent={{150,-100},{170,-80}})));
  ThermofluidStream.FlowControl.MCV mCV_cross(
    redeclare package Medium = Medium_CrossChannel,
    initM_flow=ThermofluidStream.Utilities.Types.InitializationMethods.state,
    m_flow_0=0,
    mode=ThermofluidStream.FlowControl.Internal.Types.MassflowControlValveMode.volume_flow,
    setpointFromInput=true,
    volumeFlow_set_par=0.08333) annotation (Placement(transformation(extent={{-76,-30},{-56,-10}})));
  Modelica.Blocks.Sources.Constant CrossFlowRate(k(
      quantity="VolumeFlowRate",
      unit="m3/s",
      displayUnit="m3/h") = 0.0083333333333333, y(
      quantity="VolumeFlowRate",
      unit="m3/s",
      displayUnit="m3/h")) annotation (Placement(transformation(extent={{-100,34},{-80,54}})));
  Modelica.Blocks.Sources.Constant CrossInflowTemperature(k(
      quantity="Temperature",
      unit="K",
      displayUnit="degC") = 296.35, y(
      quantity="Temperature",
      unit="K",
      displayUnit="degC")) annotation (Placement(transformation(extent={{-80,0},{-100,20}})));
  Modelica.Blocks.Sources.Constant CrossInflowPressure(k(
      quantity="AbsolutePressure",
      unit="Pa",
      displayUnit="bar") = 200000, y(
      quantity="AbsolutePressure",
      unit="Pa",
      displayUnit="bar")) annotation (Placement(transformation(extent={{-100,70},{-80,90}})));
  Modelica.Blocks.Sources.Constant CrossInflowHumidity(k(
      quantity="MassFraction",
      unit="kg/kg",
      displayUnit="g/kg") = 0.01, y(
      quantity="MassFraction",
      unit="kg/kg",
      displayUnit="g/kg")) annotation (Placement(transformation(extent={{-140,-36},{-120,-16}})));
  ThermofluidStream.Sensors.SingleSensorX X_CrossOut_Sensor(redeclare package Medium = Medium_CrossChannel, digits=5) annotation (Placement(transformation(extent={{156,-20},{180,0}})));
  ThermofluidStream.Sensors.SingleSensorSelect Enthalpy_CrossOut_Sensor(
    redeclare package Medium = Medium_CrossChannel,
    digits=0,
    quantity=ThermofluidStream.Sensors.Internal.Types.Quantities.h_Jpkg) annotation (Placement(transformation(extent={{116,-6},{142,14}})));
  ThermofluidStream.Sensors.SingleSensorSelect Temperature_CrossOut_Sensor(redeclare package Medium = Medium_CrossChannel, quantity=ThermofluidStream.Sensors.Internal.Types.Quantities.T_C)
                                                                                                                                                                                          annotation (Placement(transformation(extent={{118,-20},{138,0}})));
  ThermofluidStream.Sensors.SingleSensorSelect Pressure_CrossIn_Sensor(
    redeclare package Medium = Medium_CrossChannel,
    digits=1,
    quantity=ThermofluidStream.Sensors.Internal.Types.Quantities.p_Pa) annotation (Placement(transformation(extent={{-42,-50},{-78,-30}})));
  ThermofluidStream.Sensors.SingleSensorSelect Enthalpy_CrossIn_Sensor(
    redeclare package Medium = Medium_CrossChannel,
    digits=0,
    quantity=ThermofluidStream.Sensors.Internal.Types.Quantities.h_Jpkg) annotation (Placement(transformation(extent={{-38,-10},{-64,10}})));
  Modelica.Thermal.HeatTransfer.Celsius.PrescribedTemperature tinTemperature annotation (Placement(transformation(extent={{100,40},{80,60}})));
  Modelica.Blocks.Sources.Constant T_wall(k(
      quantity="Temperature_degC",
      unit="degC") = 5.0, y(
      quantity="Temperature_degC",
      unit="degC")) annotation (Placement(transformation(extent={{140,40},{120,60}})));
  ThermofluidStream.Boundaries.TerminalSource terminalSource[nCells](redeclare package Medium = ThermofluidStream.Media.myMedia.Water.StandardWater, each p_0=100100) annotation (Placement(transformation(extent={{-30,-50},{-10,-30}})));
  ThermofluidStream.Boundaries.Sink sink(redeclare package Medium = ThermofluidStream.Media.myMedia.Water.StandardWater, p0_par=100000) annotation (Placement(transformation(extent={{190,-70},{210,-50}})));
  ThermofluidStream.Topology.JunctionN junction_water(redeclare package Medium = ThermofluidStream.Media.myMedia.Water.StandardWater, N=nCells) annotation (Placement(transformation(extent={{80,-42},{100,-22}})));
  ThermofluidStream.Topology.JunctionN junction_air(redeclare package Medium = Medium_CrossChannel, N=nCells) annotation (Placement(transformation(extent={{80,-30},{100,-10}})));
  ThermofluidStream.Topology.SplitterN splitter_air(redeclare package Medium = Medium_CrossChannel, N=nCells) annotation (Placement(transformation(extent={{-20,-30},{0,-10}})));
  ThermofluidStream.Sensors.SingleFlowSensor MassFlow_WaterOut_Sensor(
    redeclare package Medium = Medium_Water,
    digits=3,
    quantity=ThermofluidStream.Sensors.Internal.Types.MassFlowQuantities.m_flow_gps) annotation (Placement(transformation(extent={{-10,-10},{10,10}},
        rotation=0,
        origin={126,-54})));
  ThermofluidStream.Sensors.SingleFlowSensor EnthalpyFlow_WaterOut_Sensor(
    redeclare package Medium = Medium_Water,
    digits=1,
    quantity=ThermofluidStream.Sensors.Internal.Types.MassFlowQuantities.H_flow_Jps) annotation (Placement(transformation(extent={{-10,-10},{10,10}},
        rotation=0,
        origin={166,-54})));
  ThermofluidStream.Sensors.SingleSensorSelect Temperature_WaterOut_Sensor(redeclare package Medium = Medium_Water, quantity=ThermofluidStream.Sensors.Internal.Types.Quantities.T_C) annotation (Placement(transformation(extent={{-10,-10},{10,10}},
        rotation=0,
        origin={166,-38})));
  HXutilities.Sensors.SingleSensorHumidity phi_CrossOut_Sensor(redeclare package Medium = Medium_CrossChannel, digits=2) annotation (Placement(transformation(extent={{156,-6},{176,14}})));
equation
  connect(CrossFlowRate.y, mCV_cross.setpoint_var) annotation (Line(points={{-79,44},{-66,44},{-66,-12}}, color={0,0,127}));
  connect(source_CrossChannel.outlet, mCV_cross.inlet) annotation (Line(
      points={{-82,-20},{-76,-20}},
      color={28,108,200},
      thickness=0.5));
  connect(source_CrossChannel.T0_var, CrossInflowTemperature.y) annotation (Line(points={{-94,-20},{-110,-20},{-110,10},{-101,10}},color={0,0,127}));
  connect(Pressure_CrossIn_Sensor.inlet, mCV_cross.outlet) annotation (Line(
      points={{-42,-40},{-36,-40},{-36,-20},{-56,-20}},
      color={28,108,200},
      thickness=0.5));
  connect(Enthalpy_CrossIn_Sensor.inlet, mCV_cross.outlet) annotation (Line(
      points={{-38,0},{-36,0},{-36,-20},{-56,-20}},
      color={28,108,200},
      thickness=0.5));
  connect(CrossInflowHumidity.y, source_CrossChannel.xi_var[1]) annotation (Line(points={{-119,-26},{-94,-26}}, color={0,0,127}));
  connect(tinTemperature.T, T_wall.y) annotation (Line(points={{102,50},{119,50}},
                                                                                 color={0,0,127}));
  connect(conductionElement_Cross.outlet_water, junction_water.inlets) annotation (Line(
      points={{60,-32},{80,-32}},
      color={28,108,200},
      thickness=0.5));
  connect(conductionElement_Cross.outlet, junction_air.inlets) annotation (Line(
      points={{60,-20},{80,-20}},
      color={28,108,200},
      thickness=0.5));
  connect(junction_air.outlet, Temperature_CrossOut_Sensor.inlet) annotation (Line(
      points={{100,-20},{108,-20},{108,-10},{118,-10}},
      color={28,108,200},
      thickness=0.5));
  connect(junction_air.outlet, Enthalpy_CrossOut_Sensor.inlet) annotation (Line(
      points={{100,-20},{108,-20},{108,4},{116,4}},
      color={28,108,200},
      thickness=0.5));
  connect(splitter_air.inlet, mCV_cross.outlet) annotation (Line(
      points={{-20,-20},{-56,-20}},
      color={28,108,200},
      thickness=0.5));
  connect(conductionElement_Cross.inlet, splitter_air.outlets) annotation (Line(
      points={{20,-20},{0,-20}},
      color={28,108,200},
      thickness=0.5));

  for i in 1:nCells loop
  connect(tinTemperature.port, conductionElement_Cross[i].heatPort) annotation (Line(points={{80,50},{40,50},{40,-0.4}}, color={191,0,0}));
  end for;
  connect(terminalSource.outlet, conductionElement_Cross.inlet_water) annotation (Line(
      points={{-10,-40},{4,-40},{4,-32},{20,-32}},
      color={28,108,200},
      thickness=0.5));
  connect(junction_air.outlet, X_CrossOut_Sensor.inlet) annotation (Line(
      points={{100,-20},{148,-20},{148,-10},{156,-10}},
      color={28,108,200},
      thickness=0.5));
  connect(MassFlow_WaterOut_Sensor.outlet,EnthalpyFlow_WaterOut_Sensor. inlet) annotation (Line(
      points={{136,-60},{156,-60}},
      color={28,108,200},
      thickness=0.5));
  connect(Temperature_WaterOut_Sensor.inlet,MassFlow_WaterOut_Sensor. outlet) annotation (Line(
      points={{156,-38},{144,-38},{144,-60},{136,-60}},
      color={28,108,200},
      thickness=0.5));
  connect(EnthalpyFlow_WaterOut_Sensor.outlet, sink.inlet) annotation (Line(
      points={{176,-60},{190,-60}},
      color={28,108,200},
      thickness=0.5));
  connect(junction_water.outlet, MassFlow_WaterOut_Sensor.inlet) annotation (Line(
      points={{100,-32},{108,-32},{108,-60},{116,-60}},
      color={28,108,200},
      thickness=0.5));
  connect(phi_CrossOut_Sensor.inlet, X_CrossOut_Sensor.inlet) annotation (Line(
      points={{156,4},{148,4},{148,-10},{156,-10}},
      color={28,108,200},
      thickness=0.5));
  connect(sink_CrossChannel.inlet, junction_air.outlet) annotation (Line(
      points={{190,-20},{100,-20}},
      color={28,108,200},
      thickness=0.5));
  annotation (
    Icon(coordinateSystem(preserveAspectRatio=false, extent={{-140,-100},{200,100}})),
    Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-140,-100},{200,100}}), graphics={Text(
          extent={{150,-60},{190,-82}},
          textColor={28,108,200},
          textString="Cross Flow"),
        Text(
          extent={{-3,3},{3,-3}},
          textColor={0,0,0},
          textString="g/s",
          origin={139,-53},
          rotation=0),
        Text(
          extent={{-3,3},{3,-3}},
          textColor={0,0,0},
          origin={177,3},
          rotation=0,
          textString="%%"),
        Text(
          extent={{0,3},{0,-3}},
          textColor={0,0,0},
          origin={184,-9},
          rotation=0,
          textString="kg/kg")}),
    experiment(
      StopTime=100,
      Interval=0.01,
      __Dymola_Algorithm="Esdirk23a"),
    __Dymola_Commands);
end Test_ConductionElement_Cross_CondensingAir_Array;
