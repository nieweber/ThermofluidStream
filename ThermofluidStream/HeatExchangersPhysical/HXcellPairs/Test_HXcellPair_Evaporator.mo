within ThermofluidStream.HeatExchangersPhysical.HXcellPairs;
model Test_HXcellPair_Evaporator "Test of heat exchanger cell pair: Condensing Air to Refrigerant R134a"
    extends Modelica.Icons.Example;

  replaceable package Medium_CrossChannel =
      ThermofluidStream.Media.myMedia.Air.MoistAir                        constrainedby ThermofluidStream.Media.myMedia.Interfaces.PartialMedium
                                                        annotation(choicesAllMatching = true, Dialog(group = "Medium definitions"));

  replaceable package Medium_PipeChannel = ThermofluidStream.Media.myMedia.R134a.R134a_ph
                                                                          constrainedby ThermofluidStream.Media.myMedia.Interfaces.PartialMedium
                                                        annotation(choicesAllMatching = true, Dialog(group = "Medium definitions"));

  replaceable package Medium_Water =
      ThermofluidStream.Media.myMedia.Water.StandardWater "Medium model for condensed vapor in cross channel (A)"
    annotation (choicesAllMatching=true, Dialog(group = "Medium definitions"));

  HXcellPair_CondensingAirCross_2PhasePipe hXcellPair_Evaporator(
    redeclare package Medium_CrossChannel = Medium_CrossChannel,
    redeclare package Medium_Water = Medium_Water,
    redeclare package Medium_PipeChannel = Medium_PipeChannel,
    T_0_pipe=282.15,
    T_0_cross=293.15,
    Tin(T(fixed=true))) annotation (Placement(transformation(extent={{40,-40},{80,0}})));
  ThermofluidStream.Boundaries.Source source_CrossChannel(
    redeclare package Medium = Medium_CrossChannel,
    temperatureFromInput=true,
    xiFromInput=true,
    T0_par=293.15,
    p0_par=110000) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={-82,-20})));
  ThermofluidStream.Boundaries.Source source_PipeChannel(
    redeclare package Medium = Medium_PipeChannel,
    setEnthalpy=true,
    pressureFromInput=true,
    temperatureFromInput=false,
    enthalpyFromInput=true,
    T0_par=313.15,
    p0_par=101000) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={30,-84})));
  ThermofluidStream.Boundaries.Sink sink_PipeChannel(redeclare package Medium = Medium_PipeChannel,
    pressureFromInput=true,
    p0_par=100000)                                                                                                 annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={120,86})));
  ThermofluidStream.Processes.FlowResistance flowResistance_Pipe(
    redeclare package Medium = Medium_PipeChannel,
    r=0.002,
    l=4,
    redeclare function pLoss = ThermofluidStream.Processes.Internal.FlowResistance.laminarPressureLoss) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={90,86})));
  inner ThermofluidStream.DropOfCommons dropOfCommons(m_flow_reg=1.0E-15,
                                                      assertionLevel=AssertionLevel.warning)
                                                                                annotation (Placement(transformation(extent={{130,-100},{150,-80}})));
  ThermofluidStream.FlowControl.MCV mCV_cross(
    redeclare package Medium = Medium_CrossChannel,
    initM_flow=ThermofluidStream.Utilities.Types.InitializationMethods.state,
    m_flow_0=0,
    mode=ThermofluidStream.FlowControl.Internal.Types.MassflowControlValveMode.volume_flow,
    setpointFromInput=true,
    volumeFlow_set_par=0.08333) annotation (Placement(transformation(extent={{-66,-30},{-46,-10}})));
  Modelica.Blocks.Sources.Constant CrossFlowRate(k(
      quantity="VolumeFlowRate",
      unit="m3/s",
      displayUnit="m3/h") = 0.0083333333333333, y(
      quantity="VolumeFlowRate",
      unit="m3/s",
      displayUnit="m3/h")) annotation (Placement(transformation(extent={{-90,30},{-70,50}})));
  Modelica.Blocks.Sources.Constant PipeInflowEnthalpy(k(
      quantity="SpecificEnthalpy",
      unit="J/kg",
      displayUnit="kJ/kg") = 350000, y(
      quantity="SpecificEnthalpy",
      unit="J/kg",
      displayUnit="kJ/kg")) annotation (Placement(transformation(extent={{-26,-94},{-6,-74}})));
  Modelica.Blocks.Sources.Constant PipeInflowPressure(k(
      quantity="AbsolutePressure",
      unit="Pa",
      displayUnit="bar") = 300000,  y(
      quantity="AbsolutePressure",
      unit="Pa",
      displayUnit="bar")) annotation (Placement(transformation(extent={{-26,-60},{-6,-40}})));
  Modelica.Blocks.Sources.Constant PipeMassFlowRate(k(
      quantity="MassFlowRate",
      unit="kg/s",
      displayUnit="kg/h") = 0.0055555555555556, y(
      quantity="MassFlowRate",
      unit="kg/s",
      displayUnit="kg/h")) annotation (Placement(transformation(extent={{20,54},{40,74}})));
  ThermofluidStream.FlowControl.MCV mCV_pipe(
    redeclare package Medium = Medium_PipeChannel,
    initM_flow=ThermofluidStream.Utilities.Types.InitializationMethods.state,
    m_flow_0=0,
    mode=ThermofluidStream.FlowControl.Internal.Types.MassflowControlValveMode.mass_flow,
    setpointFromInput=true,
    volumeFlow_set_par(displayUnit="l/h") = 1.9444444444444e-06,
    dp_int(start=-1000))                                         annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=90,
        origin={60,64})));
  Modelica.Blocks.Sources.Constant CrossInflowTemperature(k(
      quantity="Temperature",
      unit="K",
      displayUnit="degC") = 296.35, y(
      quantity="Temperature",
      unit="K",
      displayUnit="degC")) annotation (Placement(transformation(extent={{-130,-20},{-110,0}})));
  ThermofluidStream.Sensors.TwoPhaseSensorSelect VaporQuality_Source_Sensor(
    redeclare package Medium = Medium_PipeChannel,
    digits=3,
    quantity=ThermofluidStream.Sensors.Internal.Types.TwoPhaseQuantities.x_kgpkg) annotation (Placement(transformation(extent={{70,-82},{90,-62}})));
  ThermofluidStream.Sensors.TwoPhaseSensorSelect VaporQuality_PipeOut_Sensor(
    redeclare package Medium = Medium_PipeChannel,
    digits=3,
    quantity=ThermofluidStream.Sensors.Internal.Types.TwoPhaseQuantities.x_kgpkg) annotation (Placement(transformation(extent={{50,0},{30,20}})));
  Modelica.Blocks.Sources.Constant PipeOutflowPressure(k(
      quantity="AbsolutePressure",
      unit="Pa",
      displayUnit="bar") = 140000, y(
      quantity="AbsolutePressure",
      unit="Pa",
      displayUnit="bar")) annotation (Placement(transformation(extent={{160,76},{140,96}})));
  ThermofluidStream.Sensors.SingleSensorSelect Temperature_PipeIn_Sensor(redeclare package Medium = Medium_PipeChannel, quantity=ThermofluidStream.Sensors.Internal.Types.Quantities.T_C) annotation (Placement(transformation(extent={{70,-94},{90,-74}})));
  ThermofluidStream.Sensors.TwoPhaseSensorSelect TemperatureSaturation_PipeIn_Sensor(redeclare package Medium = Medium_PipeChannel, quantity=ThermofluidStream.Sensors.Internal.Types.TwoPhaseQuantities.T_sat_C) annotation (Placement(transformation(extent={{70,-106},{90,-86}})));
  ThermofluidStream.Sensors.SingleSensorSelect Temperature_PipeOut_Sensor(redeclare package Medium = Medium_PipeChannel, quantity=ThermofluidStream.Sensors.Internal.Types.Quantities.T_C)
                                                                                                                                                                                          annotation (Placement(transformation(extent={{70,0},{90,20}})));
  ThermofluidStream.Sensors.SingleSensorSelect Enthalpy_PipeOut_Sensor(
    redeclare package Medium = Medium_PipeChannel,
    digits=0,
    quantity=ThermofluidStream.Sensors.Internal.Types.Quantities.h_Jpkg) annotation (Placement(transformation(extent={{70,10},{96,30}})));
  ThermofluidStream.Sensors.SingleSensorSelect Pressure_PipeOut_Sensor(
    redeclare package Medium = Medium_PipeChannel,
    digits=0,
    quantity=ThermofluidStream.Sensors.Internal.Types.Quantities.p_Pa) annotation (Placement(transformation(extent={{52,10},{20,30}})));
  Modelica.Blocks.Sources.Constant CrossInflowPressure(k(
      quantity="AbsolutePressure",
      unit="Pa",
      displayUnit="bar") = 200000, y(
      quantity="AbsolutePressure",
      unit="Pa",
      displayUnit="bar")) annotation (Placement(transformation(extent={{-90,70},{-70,90}})));
  Modelica.Blocks.Sources.Ramp CrossInflowHumidity(
    height=0.03,
    duration=990,
    y(quantity="MassFraction",
      unit="kg/kg",
      displayUnit="g/kg"),
    startTime=10) annotation (Placement(transformation(extent={{-130,-94},{-110,-74}})));
  Modelica.Blocks.Sources.Constant CrossInflowHumidity1(k(
      quantity="MassFraction",
      unit="kg/kg",
      displayUnit="g/kg") = 0.015,y(
      quantity="MassFraction",
      unit="kg/kg",
      displayUnit="g/kg")) annotation (Placement(transformation(extent={{-130,-54},{-110,-34}})));
  Modelica.Blocks.Interaction.Show.RealValue MassFlowValue(
    use_numberPort=false,
    number(
      quantity="MassFlowRate",
      unit="kg/s",
      displayUnit="kg/h") = mCV_pipe.m_flow,
    significantDigits=3)                   annotation (Placement(transformation(extent={{76,64},{96,80}})));
  Modelica.Blocks.Interaction.Show.RealValue VolumeFlowValue(
    use_numberPort=false,
    number(
      quantity="VolumeFlowRate",
      unit="m3/s",
      displayUnit="l/h") = mCV_pipe.V_flow,
    significantDigits=3)                  annotation (Placement(transformation(extent={{76,52},{96,68}})));
  ThermofluidStream.Boundaries.TerminalSource terminalSource(redeclare package Medium = Medium_Water, p_0=100100)                                        annotation (Placement(transformation(extent={{2,-42},{22,-22}})));
  ThermofluidStream.Boundaries.Sink sink_CrossChannel(redeclare package Medium = Medium_CrossChannel, p0_par=101300) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={180,-20})));
  ThermofluidStream.Sensors.SingleSensorSelect Temperature_CrossOut_Sensor(redeclare package Medium = Medium_CrossChannel, quantity=ThermofluidStream.Sensors.Internal.Types.Quantities.T_C) annotation (Placement(transformation(extent={{142,-12},{162,8}})));
  ThermofluidStream.Sensors.SingleFlowSensor EnthalpyFlow_CrossOut_Sensor(
    redeclare package Medium = Medium_CrossChannel,
    digits=1,
    quantity=ThermofluidStream.Sensors.Internal.Types.MassFlowQuantities.H_flow_Jps) annotation (Placement(transformation(extent={{142,-24},{162,-4}})));
  ThermofluidStream.Sensors.SingleSensorX Humidity_CrossOut_SensorX(redeclare package Medium = Medium_CrossChannel, digits=4) annotation (Placement(transformation(extent={{138,30},{178,50}})));
  ThermofluidStream.Sensors.SingleFlowSensor MassFlow_WaterOut_Sensor(
    redeclare package Medium = Medium_Water,
    digits=3,
    quantity=ThermofluidStream.Sensors.Internal.Types.MassFlowQuantities.m_flow_gps) annotation (Placement(transformation(extent={{108,-44},{128,-24}})));
  HXutilities.Sensors.SingleSensorHumidity phi_CrossOut_Sensor(redeclare package Medium = Medium_CrossChannel, digits=2) annotation (Placement(transformation(extent={{142,16},{162,36}})));
  ThermofluidStream.Boundaries.Sink sink(redeclare package Medium = Medium_Water, p0_par=100000)                                        annotation (Placement(transformation(extent={{170,-50},{190,-30}})));
  ThermofluidStream.Sensors.SingleFlowSensor MassFlow_CrossOut_Sensor(
    redeclare package Medium = Medium_CrossChannel,
    digits=3,
    quantity=ThermofluidStream.Sensors.Internal.Types.MassFlowQuantities.m_flow_gps) annotation (Placement(transformation(extent={{100,-24},{120,-4}})));
  ThermofluidStream.Sensors.SingleFlowSensor EnthalpyFlow_WaterOut_Sensor(
    redeclare package Medium = Medium_Water,
    digits=1,
    quantity=ThermofluidStream.Sensors.Internal.Types.MassFlowQuantities.H_flow_Jps) annotation (Placement(transformation(extent={{140,-44},{160,-24}})));
  ThermofluidStream.Sensors.SingleSensorSelect Pressure_CrossIn_Sensor(
    redeclare package Medium = Medium_CrossChannel,
    digits=1,
    quantity=ThermofluidStream.Sensors.Internal.Types.Quantities.p_Pa) annotation (Placement(transformation(extent={{-26,32},{10,52}})));
  ThermofluidStream.Sensors.SingleSensorSelect Enthalpy_CrossIn_Sensor(
    redeclare package Medium = Medium_CrossChannel,
    digits=0,
    quantity=ThermofluidStream.Sensors.Internal.Types.Quantities.h_Jpkg) annotation (Placement(transformation(extent={{-26,4},{0,24}})));
  HXutilities.Sensors.SingleSensorHumidity phi_CrossIn_Sensor(redeclare package Medium = Medium_CrossChannel, digits=2) annotation (Placement(transformation(extent={{-24,18},{-2,38}})));
  ThermofluidStream.Sensors.SingleFlowSensor EnthalpyFlow_CrossIn_Sensor(
    redeclare package Medium = Medium_CrossChannel,
    digits=1,
    quantity=ThermofluidStream.Sensors.Internal.Types.MassFlowQuantities.H_flow_Jps) annotation (Placement(transformation(extent={{2,-24},{22,-4}})));
  ThermofluidStream.Sensors.SingleFlowSensor MassFlow_CrossIn_Sensor(
    redeclare package Medium = Medium_CrossChannel,
    digits=3,
    quantity=ThermofluidStream.Sensors.Internal.Types.MassFlowQuantities.m_flow_gps) annotation (Placement(transformation(extent={{-26,-24},{-6,-4}})));
  ThermofluidStream.Sensors.SingleSensorSelect Enthalpy_CrossOut_Sensor(
    redeclare package Medium = Medium_CrossChannel,
    digits=0,
    quantity=ThermofluidStream.Sensors.Internal.Types.Quantities.h_Jpkg) annotation (Placement(transformation(extent={{140,2},{166,22}})));
  ThermofluidStream.Sensors.SingleSensorSelect Temperature_CrossIn_Sensor(redeclare package Medium = Medium_CrossChannel, quantity=ThermofluidStream.Sensors.Internal.Types.Quantities.T_C) annotation (Placement(transformation(extent={{-24,-10},{-4,10}})));
  ThermofluidStream.Sensors.SingleFlowSensor EnthalpyFlow_PipeIn_Sensor(
    redeclare package Medium = Medium_PipeChannel,
    digits=1,
    quantity=ThermofluidStream.Sensors.Internal.Types.MassFlowQuantities.H_flow_Jps) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=90,
        origin={54,-60})));
  ThermofluidStream.Sensors.SingleFlowSensor EnthalpyFlow_PipeOut_Sensor(
    redeclare package Medium = Medium_PipeChannel,
    digits=1,
    quantity=ThermofluidStream.Sensors.Internal.Types.MassFlowQuantities.H_flow_Jps) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=90,
        origin={54,38})));
  ThermofluidStream.Sensors.SingleSensorSelect Temperature_WaterOut_Sensor(redeclare package Medium = Medium_Water, quantity=ThermofluidStream.Sensors.Internal.Types.Quantities.T_C) annotation (Placement(transformation(extent={{142,-62},{162,-42}})));
equation
  connect(flowResistance_Pipe.outlet, sink_PipeChannel.inlet) annotation (Line(
      points={{100,86},{110,86}},
      color={28,108,200},
      thickness=0.5));
  connect(CrossFlowRate.y, mCV_cross.setpoint_var) annotation (Line(points={{-69,40},{-56,40},{-56,-12}}, color={0,0,127}));
  connect(source_CrossChannel.outlet, mCV_cross.inlet) annotation (Line(
      points={{-72,-20},{-66,-20}},
      color={28,108,200},
      thickness=0.5));
  connect(source_PipeChannel.p0_var, PipeInflowPressure.y) annotation (Line(points={{28,-78},{28,-50},{-5,-50}},              color={0,0,127}));
  connect(PipeMassFlowRate.y,mCV_pipe. setpoint_var) annotation (Line(points={{41,64},{52,64}},  color={0,0,127}));
  connect(mCV_pipe.outlet, flowResistance_Pipe.inlet) annotation (Line(
      points={{60,74},{60,86},{80,86}},
      color={28,108,200},
      thickness=0.5));
  connect(source_CrossChannel.T0_var, CrossInflowTemperature.y) annotation (Line(points={{-84,-20},{-100,-20},{-100,-10},{-109,-10}},
                                                                                                                                   color={0,0,127}));
  connect(PipeInflowEnthalpy.y, source_PipeChannel.h0_var) annotation (Line(points={{-5,-84},{28,-84}},                       color={0,0,127}));
  connect(VaporQuality_Source_Sensor.inlet, source_PipeChannel.outlet) annotation (Line(
      points={{70,-72},{60,-72},{60,-84},{40,-84}},
      color={28,108,200},
      thickness=0.5));
  connect(hXcellPair_Evaporator.outlet_Pipe, VaporQuality_PipeOut_Sensor.inlet) annotation (Line(
      points={{60,0},{60,10},{50,10}},
      color={28,108,200},
      thickness=0.5));
  connect(sink_PipeChannel.p0_var, PipeOutflowPressure.y) annotation (Line(points={{122,86},{139,86}},
                                                                                                     color={0,0,127}));
  connect(TemperatureSaturation_PipeIn_Sensor.inlet, source_PipeChannel.outlet) annotation (Line(
      points={{70,-96},{60,-96},{60,-84},{40,-84}},
      color={28,108,200},
      thickness=0.5));
  connect(Temperature_PipeIn_Sensor.inlet, source_PipeChannel.outlet) annotation (Line(
      points={{70,-84},{40,-84}},
      color={28,108,200},
      thickness=0.5));
  connect(hXcellPair_Evaporator.outlet_Pipe, Temperature_PipeOut_Sensor.inlet) annotation (Line(
      points={{60,0},{60,10},{70,10}},
      color={28,108,200},
      thickness=0.5));
  connect(Enthalpy_PipeOut_Sensor.inlet, hXcellPair_Evaporator.outlet_Pipe) annotation (Line(
      points={{70,20},{60,20},{60,0}},
      color={28,108,200},
      thickness=0.5));
  connect(hXcellPair_Evaporator.outlet_Pipe, Pressure_PipeOut_Sensor.inlet) annotation (Line(
      points={{60,0},{60,20},{52,20}},
      color={28,108,200},
      thickness=0.5));
  connect(CrossInflowHumidity1.y, source_CrossChannel.xi_var[1]) annotation (Line(points={{-109,-44},{-100,-44},{-100,-26},{-84,-26}},
                                                                                                                 color={0,0,127}));
  connect(terminalSource.outlet, hXcellPair_Evaporator.inlet_cross_water) annotation (Line(
      points={{22,-32},{32,-32},{32,-28},{40,-28}},
      color={28,108,200},
      thickness=0.5));
  connect(EnthalpyFlow_CrossOut_Sensor.outlet,sink_CrossChannel. inlet) annotation (Line(
      points={{162,-20},{170,-20}},
      color={28,108,200},
      thickness=0.5));
  connect(MassFlow_CrossOut_Sensor.outlet, EnthalpyFlow_CrossOut_Sensor.inlet) annotation (Line(
      points={{120,-20},{142,-20}},
      color={28,108,200},
      thickness=0.5));
  connect(MassFlow_CrossOut_Sensor.outlet, Temperature_CrossOut_Sensor.inlet) annotation (Line(
      points={{120,-20},{130,-20},{130,-2},{142,-2}},
      color={28,108,200},
      thickness=0.5));
  connect(phi_CrossOut_Sensor.inlet, MassFlow_CrossOut_Sensor.outlet) annotation (Line(
      points={{142,26},{130,26},{130,-20},{120,-20}},
      color={28,108,200},
      thickness=0.5));
  connect(Humidity_CrossOut_SensorX.inlet, MassFlow_CrossOut_Sensor.outlet) annotation (Line(
      points={{138,40},{130,40},{130,-20},{120,-20}},
      color={28,108,200},
      thickness=0.5));
  connect(EnthalpyFlow_WaterOut_Sensor.outlet, sink.inlet) annotation (Line(
      points={{160,-40},{170,-40}},
      color={28,108,200},
      thickness=0.5));
  connect(MassFlow_WaterOut_Sensor.outlet, EnthalpyFlow_WaterOut_Sensor.inlet) annotation (Line(
      points={{128,-40},{140,-40}},
      color={28,108,200},
      thickness=0.5));
  connect(hXcellPair_Evaporator.outlet_cross_water, MassFlow_WaterOut_Sensor.inlet) annotation (Line(
      points={{80,-28},{94,-28},{94,-40},{108,-40}},
      color={28,108,200},
      thickness=0.5));
  connect(hXcellPair_Evaporator.outlet_Cross, MassFlow_CrossOut_Sensor.inlet) annotation (Line(
      points={{80,-20},{100,-20}},
      color={28,108,200},
      thickness=0.5));
  connect(MassFlow_CrossIn_Sensor.outlet, EnthalpyFlow_CrossIn_Sensor.inlet) annotation (Line(
      points={{-6,-20},{2,-20}},
      color={28,108,200},
      thickness=0.5));
  connect(EnthalpyFlow_CrossIn_Sensor.outlet, hXcellPair_Evaporator.inlet_Cross) annotation (Line(
      points={{22,-20},{40,-20}},
      color={28,108,200},
      thickness=0.5));
  connect(mCV_cross.outlet, MassFlow_CrossIn_Sensor.inlet) annotation (Line(
      points={{-46,-20},{-26,-20}},
      color={28,108,200},
      thickness=0.5));
  connect(Enthalpy_CrossIn_Sensor.inlet, mCV_cross.outlet) annotation (Line(
      points={{-26,14},{-36,14},{-36,-20},{-46,-20}},
      color={28,108,200},
      thickness=0.5));
  connect(phi_CrossIn_Sensor.inlet, mCV_cross.outlet) annotation (Line(
      points={{-24,28},{-36,28},{-36,-20},{-46,-20}},
      color={28,108,200},
      thickness=0.5));
  connect(Pressure_CrossIn_Sensor.inlet, mCV_cross.outlet) annotation (Line(
      points={{-26,42},{-36,42},{-36,-20},{-46,-20}},
      color={28,108,200},
      thickness=0.5));
  connect(Enthalpy_CrossOut_Sensor.inlet, MassFlow_CrossOut_Sensor.outlet) annotation (Line(
      points={{140,12},{130,12},{130,-20},{120,-20}},
      color={28,108,200},
      thickness=0.5));
  connect(Temperature_CrossIn_Sensor.inlet, mCV_cross.outlet) annotation (Line(
      points={{-24,0},{-36,0},{-36,-20},{-46,-20}},
      color={28,108,200},
      thickness=0.5));
  connect(hXcellPair_Evaporator.inlet_Pipe, EnthalpyFlow_PipeIn_Sensor.outlet) annotation (Line(
      points={{60,-40},{60,-50}},
      color={28,108,200},
      thickness=0.5));
  connect(EnthalpyFlow_PipeIn_Sensor.inlet, source_PipeChannel.outlet) annotation (Line(
      points={{60,-70},{60,-84},{40,-84}},
      color={28,108,200},
      thickness=0.5));
  connect(mCV_pipe.inlet, EnthalpyFlow_PipeOut_Sensor.outlet) annotation (Line(
      points={{60,54},{60,48}},
      color={28,108,200},
      thickness=0.5));
  connect(EnthalpyFlow_PipeOut_Sensor.inlet, hXcellPair_Evaporator.outlet_Pipe) annotation (Line(
      points={{60,28},{60,0}},
      color={28,108,200},
      thickness=0.5));
  connect(Temperature_WaterOut_Sensor.inlet, MassFlow_WaterOut_Sensor.outlet) annotation (Line(
      points={{142,-52},{134,-52},{134,-40},{128,-40}},
      color={28,108,200},
      thickness=0.5));
  annotation (
    Icon(coordinateSystem(preserveAspectRatio=false, extent={{-140,-100},{180,100}})),
    Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-140,-100},{180,100}}), graphics={
        Text(
          extent={{116,72},{156,50}},
          textColor={0,0,0},
          textString="Pipe Flow"),
        Text(
          extent={{-66,-90},{-34,-96}},
          textColor={28,108,200},
          textString="2-phase region:
10 bar: 250 .. 420 kJ/kg"),
        Text(
          extent={{18,92},{46,80}},
          textColor={28,108,200},
          textString="GUNT: 10.5 kg/s"),
        Text(
          extent={{136,-60},{176,-80}},
          textColor={0,0,0},
          textString="Cross Flow"),
        Text(
          extent={{128,-30},{134,-36}},
          textColor={0,0,0},
          textString="g/s"),
        Text(
          extent={{-8,-10},{-2,-16}},
          textColor={0,0,0},
          textString="g/s"),
        Text(
          extent={{118,-10},{124,-16}},
          textColor={0,0,0},
          textString="g/s")}),
    experiment(
      StopTime=100,
      Interval=0.01,
      __Dymola_Algorithm="Esdirk23a"),
    __Dymola_Commands);
end Test_HXcellPair_Evaporator;
