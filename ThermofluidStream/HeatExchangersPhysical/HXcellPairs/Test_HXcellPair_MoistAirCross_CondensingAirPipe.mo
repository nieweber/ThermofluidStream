within ThermofluidStream.HeatExchangersPhysical.HXcellPairs;
model Test_HXcellPair_MoistAirCross_CondensingAirPipe
    extends Modelica.Icons.Example;

  replaceable package Medium_CrossChannel =
      ThermofluidStream.Media.myMedia.Air.MoistAir                        constrainedby ThermofluidStream.Media.myMedia.Interfaces.PartialMedium
                                                        annotation(choicesAllMatching = true);

  package Medium_PipeChannel = ThermofluidStream.Media.myMedia.Air.MoistAir;

//  replaceable package Medium_PipeChannel = ThermofluidStream.Media.myMedia.Air.MoistAir
//                                                                          constrainedby ThermofluidStream.Media.myMedia.Interfaces.PartialMedium
//                                                        annotation(choicesAllMatching = true);

  replaceable package Medium_Water =
      ThermofluidStream.Media.myMedia.Water.StandardWater "Medium model for condensed vapor in pipe channel"
    annotation (choicesAllMatching=true, Dialog(group = "Medium definitions"));

  HXcellPairs.HXcellPair_1PhaseCross_CondensingAirPipe hXcellPair(
    redeclare package Medium_CrossChannel = Medium_CrossChannel,
    T_pipe_initial=PipeInflowTemperature.k,
    T_cross_initial=CrossInflowTemperature.k) annotation (Placement(transformation(extent={{-20,-40},{20,0}})));

  ThermofluidStream.Boundaries.Source source_CrossChannel(
    redeclare package Medium = Medium_CrossChannel,
    temperatureFromInput=true,
    xiFromInput=true,
    T0_par=293.15,
    p0_par=100010) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={-110,-20})));
  ThermofluidStream.Boundaries.Source source_PipeChannel(
    redeclare package Medium = Medium_PipeChannel,
    pressureFromInput=true,
    temperatureFromInput=true,
    xiFromInput=true,
    T0_par=313.15,
    p0_par=101000) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={-98,-60})));
  ThermofluidStream.Boundaries.Sink sink_CrossChannel(redeclare package Medium = Medium_CrossChannel, p0_par=100000) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={154,-20})));
  ThermofluidStream.Boundaries.Sink sink_PipeChannel(redeclare package Medium = Medium_PipeChannel, p0_par=100000) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={110,80})));
  ThermofluidStream.Processes.FlowResistance flowResistance_Cross(
    redeclare package Medium = Medium_CrossChannel,
    r=0.07,
    l=0.2,
    redeclare function pLoss = ThermofluidStream.Processes.Internal.FlowResistance.laminarTurbulentPressureLossHaaland) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={122,-20})));
  inner ThermofluidStream.DropOfCommons dropOfCommons(m_flow_reg=1E-15,
                                                      assertionLevel=AssertionLevel.warning)
                                                                                annotation (Placement(transformation(extent={{60,-100},{80,-80}})));
  ThermofluidStream.FlowControl.MCV mCV_cross(
    redeclare package Medium = Medium_CrossChannel,
    mode=ThermofluidStream.FlowControl.Internal.Types.MassflowControlValveMode.volume_flow,
    setpointFromInput=true,
    volumeFlow_set_par=0.08333) annotation (Placement(transformation(extent={{-90,-30},{-70,-10}})));
  Modelica.Blocks.Sources.Constant CrossFlowRate(k(
      quantity="VolumeFlowRate",
      unit="m3/s",
      displayUnit="m3/h") = 0.0083333333333333, y(
      quantity="VolumeFlowRate",
      unit="m3/s",
      displayUnit="m3/h")) annotation (Placement(transformation(extent={{-120,0},{-100,20}})));
  Modelica.Blocks.Sources.Ramp     PipeInflowHumidity(
    height=0.03,
    duration=990,                   y(
      quantity="MassFraction",
      unit="kg/kg",
      displayUnit="g/kg"),
    startTime=10)          annotation (Placement(transformation(extent={{-180,-90},{-160,-70}})));
  Modelica.Blocks.Sources.Constant PipeInflowPressure(k(
      quantity="AbsolutePressure",
      unit="Pa",
      displayUnit="bar") = 200000,  y(
      quantity="AbsolutePressure",
      unit="Pa",
      displayUnit="bar")) annotation (Placement(transformation(extent={{-180,-60},{-160,-40}})));
  Modelica.Blocks.Sources.Constant PipeMassFlowRate(k(
      quantity="MassFlowRate",
      unit="kg/s",
      displayUnit="kg/h") = 0.0011111111111111, y(
      quantity="MassFlowRate",
      unit="kg/s",
      displayUnit="kg/h")) annotation (Placement(transformation(extent={{10,100},{30,120}})));
  ThermofluidStream.FlowControl.MCV mCV_pipe(
    redeclare package Medium = Medium_PipeChannel,
    mode=ThermofluidStream.FlowControl.Internal.Types.MassflowControlValveMode.mass_flow,
    setpointFromInput=true,
    volumeFlow_set_par(displayUnit="l/h") = 1.9444444444444e-06,
    dp_int(start=-1000))                                         annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={56,80})));
  Modelica.Blocks.Sources.Constant CrossInflowTemperature(k(
      quantity="Temperature",
      unit="K",
      displayUnit="degC") = 288.15, y(
      quantity="Temperature",
      unit="K",
      displayUnit="degC")) annotation (Placement(transformation(extent={{-180,0},{-160,20}})));
  Modelica.Blocks.Sources.Constant PipeInflowTemperature(k(
      quantity="Temperature",
      unit="K",
      displayUnit="degC") = 323.15, y(
      quantity="Temperature",
      unit="K",
      displayUnit="degC")) annotation (Placement(transformation(extent={{-140,-70},{-120,-50}})));
  ThermofluidStream.Sensors.SingleSensorSelect Temperature_PipeOut_Sensor(redeclare package Medium = Medium_PipeChannel, quantity=ThermofluidStream.Sensors.Internal.Types.Quantities.T_C)
                                                                                                                                                                                          annotation (Placement(transformation(extent={{16,8},{36,28}})));
  ThermofluidStream.Sensors.SingleSensorSelect Enthalpy_PipeOut_Sensor(
    redeclare package Medium = Medium_PipeChannel,
    digits=0,
    quantity=ThermofluidStream.Sensors.Internal.Types.Quantities.h_Jpkg) annotation (Placement(transformation(extent={{16,18},{42,38}})));
  ThermofluidStream.Sensors.SingleSensorSelect Enthalpy_PipeIn_Sensor(
    redeclare package Medium = Medium_PipeChannel,
    digits=0,
    quantity=ThermofluidStream.Sensors.Internal.Types.Quantities.h_Jpkg) annotation (Placement(transformation(extent={{10,-84},{36,-64}})));
  Modelica.Blocks.Interaction.Show.RealValue MassFlowValue(
    use_numberPort=false,
    number(
      quantity="MassFlowRate",
      unit="kg/s",
      displayUnit="kg/h") = mCV_pipe.m_flow,
    significantDigits=3)                   annotation (Placement(transformation(extent={{46,60},{66,76}})));
  Modelica.Blocks.Interaction.Show.RealValue VolumeFlowValue(
    use_numberPort=false,
    number(
      quantity="VolumeFlowRate",
      unit="m3/s",
      displayUnit="l/h") = mCV_pipe.V_flow,
    significantDigits=3)                  annotation (Placement(transformation(extent={{46,48},{66,64}})));
  Modelica.Blocks.Sources.Constant PipeInflowHumidity1(k(
      quantity="MassFraction",
      unit="kg/kg",
      displayUnit="g/kg") = 0.01, y(
      quantity="MassFraction",
      unit="kg/kg",
      displayUnit="g/kg")) annotation (Placement(transformation(extent={{-218,-90},{-198,-70}})));
  ThermofluidStream.Boundaries.TerminalSource terminalSource_water(redeclare package Medium = Medium_Water, p_0=100100) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=90,
        origin={16,-64})));
  ThermofluidStream.Sensors.SingleFlowSensor MassFlow_WaterOut_Sensor(
    redeclare package Medium = Medium_Water,
    digits=3,
    quantity=ThermofluidStream.Sensors.Internal.Types.MassFlowQuantities.m_flow_gps) annotation (Placement(transformation(extent={{-10,-10},{10,10}},
        rotation=0,
        origin={60,16})));
  ThermofluidStream.Boundaries.Sink sink_water(redeclare package Medium = Medium_Water, p0_par=100000) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={136,10})));
  ThermofluidStream.Sensors.SingleFlowSensor EnthalpyFlow_WaterOut_Sensor(
    redeclare package Medium = Medium_Water,
    digits=1,
    quantity=ThermofluidStream.Sensors.Internal.Types.MassFlowQuantities.H_flow_Jps) annotation (Placement(transformation(extent={{-10,-10},{10,10}},
        rotation=0,
        origin={100,16})));
  ThermofluidStream.Sensors.SingleSensorSelect Temperature_WaterOut_Sensor(redeclare package Medium = Medium_Water, quantity=ThermofluidStream.Sensors.Internal.Types.Quantities.T_C) annotation (Placement(transformation(extent={{-10,-10},{10,10}},
        rotation=0,
        origin={100,32})));
  ThermofluidStream.Sensors.SingleFlowSensor EnthalpyFlow_PipeIn_Sensor(
    redeclare package Medium = Medium_PipeChannel,
    digits=1,
    quantity=ThermofluidStream.Sensors.Internal.Types.MassFlowQuantities.H_flow_Jps) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={-30,-54})));
  ThermofluidStream.Sensors.SingleFlowSensor EnthalpyFlow_PipeOut_Sensor(
    redeclare package Medium = Medium_PipeChannel,
    digits=1,
    quantity=ThermofluidStream.Sensors.Internal.Types.MassFlowQuantities.H_flow_Jps) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={20,86})));
  ThermofluidStream.Sensors.SingleFlowSensor EnthalpyFlow_CrossIn_Sensor(
    redeclare package Medium = Medium_CrossChannel,
    digits=1,
    quantity=ThermofluidStream.Sensors.Internal.Types.MassFlowQuantities.H_flow_Jps) annotation (Placement(transformation(extent={{-60,-24},{-40,-4}})));
  ThermofluidStream.Sensors.SingleFlowSensor EnthalpyFlow_CrossOut_Sensor(
    redeclare package Medium = Medium_CrossChannel,
    digits=1,
    quantity=ThermofluidStream.Sensors.Internal.Types.MassFlowQuantities.H_flow_Jps) annotation (Placement(transformation(extent={{80,-24},{100,-4}})));
  ThermofluidStream.Sensors.SingleSensorSelect Temperature_CrossOut_Sensor(redeclare package Medium = Medium_CrossChannel, quantity=ThermofluidStream.Sensors.Internal.Types.Quantities.T_C) annotation (Placement(transformation(extent={{80,-46},{100,-26}})));
  HXutilities.Sensors.SingleSensorHumidity phi_CrossOut_Sensor(redeclare package Medium = Medium_CrossChannel, digits=2) annotation (Placement(transformation(extent={{80,-58},{100,-38}})));
  ThermofluidStream.Sensors.SingleFlowSensor MassFlow_CrossOut_Sensor(
    redeclare package Medium = Medium_CrossChannel,
    digits=3,
    quantity=ThermofluidStream.Sensors.Internal.Types.MassFlowQuantities.m_flow_gps) annotation (Placement(transformation(extent={{40,-24},{60,-4}})));
  Modelica.Blocks.Sources.Constant CrossInflowHumidity(k(
      quantity="MassFraction",
      unit="kg/kg",
      displayUnit="g/kg") = 0.01, y(
      quantity="MassFraction",
      unit="kg/kg",
      displayUnit="g/kg")) annotation (Placement(transformation(extent={{-180,-30},{-160,-10}})));
  HXutilities.Sensors.SingleSensorHumidity phi_CrossIn_Sensor(redeclare package Medium = Medium_CrossChannel, digits=2) annotation (Placement(transformation(extent={{-60,-12},{-38,8}})));
  ThermofluidStream.Sensors.SingleFlowSensor MassFlow_PipeIn_Sensor(
    redeclare package Medium = Medium_PipeChannel,
    digits=3,
    quantity=ThermofluidStream.Sensors.Internal.Types.MassFlowQuantities.m_flow_gps) annotation (Placement(transformation(extent={{-70,-64},{-50,-44}})));
  ThermofluidStream.Sensors.SingleFlowSensor MassFlow_PipeOut_Sensor(
    redeclare package Medium = Medium_PipeChannel,
    digits=3,
    quantity=ThermofluidStream.Sensors.Internal.Types.MassFlowQuantities.m_flow_gps) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=90,
        origin={-6,50})));
  ThermofluidStream.Sensors.SingleSensorSelect Temperature_PipeIn_Sensor(redeclare package Medium = Medium_PipeChannel, quantity=ThermofluidStream.Sensors.Internal.Types.Quantities.T_C) annotation (Placement(transformation(extent={{12,-96},{32,-76}})));
  ThermofluidStream.Sensors.SingleSensorX Humidity_PipeOut_SensorX(redeclare package Medium = Medium_CrossChannel, digits=4) annotation (Placement(transformation(extent={{-12,8},{-52,28}})));
  ThermofluidStream.Sensors.SingleSensorX Humidity_PipeIn_SensorX(redeclare package Medium = Medium_CrossChannel, digits=4) annotation (Placement(transformation(extent={{-10,-96},{-50,-76}})));
  HXutilities.Sensors.SingleSensorHumidity phi_PipeIn_Sensor(redeclare package Medium = Medium_PipeChannel, digits=2) annotation (Placement(transformation(extent={{-14,-84},{-36,-64}})));
  HXutilities.Sensors.SingleSensorHumidity phi_PipeOut_Sensor(redeclare package Medium = Medium_PipeChannel, digits=2) annotation (Placement(transformation(extent={{-16,20},{-36,40}})));
equation
  connect(flowResistance_Cross.outlet, sink_CrossChannel.inlet) annotation (Line(
      points={{132,-20},{144,-20}},
      color={28,108,200},
      thickness=0.5));
  connect(CrossFlowRate.y, mCV_cross.setpoint_var) annotation (Line(points={{-99,10},{-80,10},{-80,-12}}, color={0,0,127}));
  connect(source_CrossChannel.outlet, mCV_cross.inlet) annotation (Line(
      points={{-100,-20},{-90,-20}},
      color={28,108,200},
      thickness=0.5));
  connect(source_PipeChannel.p0_var, PipeInflowPressure.y) annotation (Line(points={{-100,-54},{-110,-54},{-110,-40},{-144,-40},{-144,-50},{-159,-50}},
                                                                                                                              color={0,0,127}));
  connect(PipeMassFlowRate.y,mCV_pipe. setpoint_var) annotation (Line(points={{31,110},{56,110},{56,88}},
                                                                                                 color={0,0,127}));
  connect(source_CrossChannel.T0_var, CrossInflowTemperature.y) annotation (Line(points={{-112,-20},{-140,-20},{-140,10},{-159,10}},
                                                                                                                                   color={0,0,127}));
  connect(PipeInflowHumidity.y, source_PipeChannel.xi_var[1]) annotation (Line(points={{-159,-80},{-110,-80},{-110,-66},{-100,-66}},
                                                                                                                                   color={0,0,127}));
  connect(PipeInflowTemperature.y, source_PipeChannel.T0_var) annotation (Line(points={{-119,-60},{-100,-60}},                   color={0,0,127}));
  connect(hXcellPair.outlet_Pipe, Temperature_PipeOut_Sensor.inlet) annotation (Line(
      points={{0,0},{0,18},{16,18}},
      color={28,108,200},
      thickness=0.5));
  connect(hXcellPair.outlet_Pipe, Enthalpy_PipeOut_Sensor.inlet) annotation (Line(
      points={{0,0},{0,18},{10,18},{10,28},{16,28}},
      color={28,108,200},
      thickness=0.5));
  connect(terminalSource_water.outlet, hXcellPair.inlet_pipe_water) annotation (Line(
      points={{16,-54},{16,-48},{8,-48},{8,-40}},
      color={28,108,200},
      thickness=0.5));
  connect(EnthalpyFlow_WaterOut_Sensor.outlet, sink_water.inlet) annotation (Line(
      points={{110,10},{126,10}},
      color={28,108,200},
      thickness=0.5));
  connect(MassFlow_WaterOut_Sensor.outlet,EnthalpyFlow_WaterOut_Sensor. inlet) annotation (Line(
      points={{70,10},{90,10}},
      color={28,108,200},
      thickness=0.5));
  connect(Temperature_WaterOut_Sensor.inlet,MassFlow_WaterOut_Sensor. outlet) annotation (Line(
      points={{90,32},{78,32},{78,10},{70,10}},
      color={28,108,200},
      thickness=0.5));
  connect(hXcellPair.outlet_pipe_water, MassFlow_WaterOut_Sensor.inlet) annotation (Line(
      points={{8,-0.4},{8,10},{50,10}},
      color={28,108,200},
      thickness=0.5));
  connect(EnthalpyFlow_PipeIn_Sensor.outlet, hXcellPair.inlet_Pipe) annotation (Line(
      points={{-20,-60},{0,-60},{0,-40}},
      color={28,108,200},
      thickness=0.5));
  connect(Enthalpy_PipeIn_Sensor.inlet, EnthalpyFlow_PipeIn_Sensor.outlet) annotation (Line(
      points={{10,-74},{0,-74},{0,-60},{-20,-60}},
      color={28,108,200},
      thickness=0.5));
  connect(mCV_cross.outlet, EnthalpyFlow_CrossIn_Sensor.inlet) annotation (Line(
      points={{-70,-20},{-60,-20}},
      color={28,108,200},
      thickness=0.5));
  connect(EnthalpyFlow_CrossIn_Sensor.outlet, hXcellPair.inlet_Cross) annotation (Line(
      points={{-40,-20},{-20,-20}},
      color={28,108,200},
      thickness=0.5));
  connect(EnthalpyFlow_CrossOut_Sensor.outlet, flowResistance_Cross.inlet) annotation (Line(
      points={{100,-20},{112,-20}},
      color={28,108,200},
      thickness=0.5));
  connect(hXcellPair.outlet_Cross, MassFlow_CrossOut_Sensor.inlet) annotation (Line(
      points={{20,-20},{40,-20}},
      color={28,108,200},
      thickness=0.5));
  connect(MassFlow_CrossOut_Sensor.outlet, EnthalpyFlow_CrossOut_Sensor.inlet) annotation (Line(
      points={{60,-20},{80,-20}},
      color={28,108,200},
      thickness=0.5));
  connect(Temperature_CrossOut_Sensor.inlet, MassFlow_CrossOut_Sensor.outlet) annotation (Line(
      points={{80,-36},{70,-36},{70,-20},{60,-20}},
      color={28,108,200},
      thickness=0.5));
  connect(phi_CrossOut_Sensor.inlet, MassFlow_CrossOut_Sensor.outlet) annotation (Line(
      points={{80,-48},{70,-48},{70,-20},{60,-20}},
      color={28,108,200},
      thickness=0.5));
  connect(CrossInflowHumidity.y, source_CrossChannel.xi_var[1]) annotation (Line(points={{-159,-20},{-148,-20},{-148,-26},{-112,-26}}, color={0,0,127}));
  connect(phi_CrossIn_Sensor.inlet, mCV_cross.outlet) annotation (Line(
      points={{-60,-2},{-66,-2},{-66,-20},{-70,-20}},
      color={28,108,200},
      thickness=0.5));
  connect(MassFlow_PipeIn_Sensor.outlet, EnthalpyFlow_PipeIn_Sensor.inlet) annotation (Line(
      points={{-50,-60},{-40,-60}},
      color={28,108,200},
      thickness=0.5));
  connect(source_PipeChannel.outlet, MassFlow_PipeIn_Sensor.inlet) annotation (Line(
      points={{-88,-60},{-70,-60}},
      color={28,108,200},
      thickness=0.5));
  connect(MassFlow_PipeOut_Sensor.inlet, hXcellPair.outlet_Pipe) annotation (Line(
      points={{0,40},{0,0}},
      color={28,108,200},
      thickness=0.5));
  connect(EnthalpyFlow_PipeOut_Sensor.inlet, MassFlow_PipeOut_Sensor.outlet) annotation (Line(
      points={{10,80},{0,80},{0,60}},
      color={28,108,200},
      thickness=0.5));
  connect(EnthalpyFlow_PipeOut_Sensor.outlet, mCV_pipe.inlet) annotation (Line(
      points={{30,80},{46,80}},
      color={28,108,200},
      thickness=0.5));
  connect(Temperature_PipeIn_Sensor.inlet, EnthalpyFlow_PipeIn_Sensor.outlet) annotation (Line(
      points={{12,-86},{0,-86},{0,-60},{-20,-60}},
      color={28,108,200},
      thickness=0.5));
  connect(Humidity_PipeOut_SensorX.inlet, hXcellPair.outlet_Pipe) annotation (Line(
      points={{-12,18},{0,18},{0,0}},
      color={28,108,200},
      thickness=0.5));
  connect(Humidity_PipeIn_SensorX.inlet, EnthalpyFlow_PipeIn_Sensor.outlet) annotation (Line(
      points={{-10,-86},{0,-86},{0,-60},{-20,-60}},
      color={28,108,200},
      thickness=0.5));
  connect(phi_PipeIn_Sensor.inlet, EnthalpyFlow_PipeIn_Sensor.outlet) annotation (Line(
      points={{-14,-74},{0,-74},{0,-60},{-20,-60}},
      color={28,108,200},
      thickness=0.5));
  connect(phi_PipeOut_Sensor.inlet, hXcellPair.outlet_Pipe) annotation (Line(
      points={{-16,30},{0,30},{0,0}},
      color={28,108,200},
      thickness=0.5));
  connect(mCV_pipe.outlet, sink_PipeChannel.inlet) annotation (Line(
      points={{66,80},{100,80}},
      color={28,108,200},
      thickness=0.5));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-180,-100},{140,120}})),
                                                                 Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-180,-100},{140,120}}),
                                                                                                                      graphics={Text(
          extent={{-52,94},{-12,72}},
          textColor={0,0,0},
          textString="Pipe Flow"), Text(
          extent={{92,-52},{132,-74}},
          textColor={0,0,0},
          textString="Cross Flow"),
        Text(
          extent={{-3,3},{3,-3}},
          textColor={0,0,0},
          textString="g/s",
          origin={73,17},
          rotation=0),
        Text(
          extent={{-136,-4},{-112,-18}},
          textColor={0,0,0},
          textString="Moist Air"),
        Text(
          extent={{118,32},{142,16}},
          textColor={0,0,0},
          textString="Liquid Water"),
        Text(
          extent={{-102,-70},{-78,-84}},
          textColor={0,0,0},
          textString="Moist Air"),
        Text(
          extent={{144,-30},{160,-42}},
          textColor={0,0,0},
          textString="Moist Air"),
        Text(
          extent={{-3,3},{3,-3}},
          textColor={0,0,0},
          textString="g/s",
          origin={63,-13},
          rotation=0),
        Text(
          extent={{-3,3},{3,-3}},
          textColor={0,0,0},
          origin={101,-49},
          rotation=0,
          textString="%%"),
        Text(
          extent={{-3,3},{3,-3}},
          textColor={0,0,0},
          origin={-37,-3},
          rotation=0,
          textString="%%"),
        Text(
          extent={{-3,3},{3,-3}},
          textColor={0,0,0},
          textString="g/s",
          origin={-49,-53},
          rotation=0),
        Text(
          extent={{-3,3},{3,-3}},
          textColor={0,0,0},
          textString="g/s",
          origin={-7,61},
          rotation=0),
        Text(
          extent={{-3,3},{3,-3}},
          textColor={0,0,0},
          origin={-35,-75},
          rotation=0,
          textString="%%"),
        Text(
          extent={{-3,3},{3,-3}},
          textColor={0,0,0},
          origin={-39,29},
          rotation=0,
          textString="%%"),
        Text(
          extent={{100,68},{116,56}},
          textColor={0,0,0},
          textString="Moist Air"),
        Text(
          extent={{0,3},{0,-3}},
          textColor={0,0,0},
          origin={-52,-87},
          rotation=0,
          textString="kg/kg"),
        Text(
          extent={{0,3},{0,-3}},
          textColor={0,0,0},
          origin={-52,19},
          rotation=0,
          textString="kg/kg")}),
    experiment(
      StopTime=1000,
      Interval=0.01,
      __Dymola_Algorithm="Dassl"),
    __Dymola_Commands(file="Scripts/Plots_HXcellPair.mos" "Plots_HXcellPair"));
end Test_HXcellPair_MoistAirCross_CondensingAirPipe;
