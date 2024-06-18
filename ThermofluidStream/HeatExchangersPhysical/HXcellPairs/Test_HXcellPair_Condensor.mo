within ThermofluidStream.HeatExchangersPhysical.HXcellPairs;
model Test_HXcellPair_Condensor "Test of heat exchanger cell pair: moist air to Refrigerant R134a"
    extends Modelica.Icons.Example;

  replaceable package Medium_CrossChannel =
      ThermofluidStream.Media.myMedia.Air.MoistAir                        constrainedby ThermofluidStream.Media.myMedia.Interfaces.PartialMedium
                                                        annotation(choicesAllMatching = true);

  replaceable package Medium_PipeChannel = ThermofluidStream.Media.myMedia.R134a.R134a_ph
                                                                          constrainedby ThermofluidStream.Media.myMedia.Interfaces.PartialMedium
                                                        annotation(choicesAllMatching = true);

  HXcellPairs.HXcellPair_1PhaseCross_2PhasePipe hXcellPair_Refrig(
    redeclare package Medium_CrossChannel = Medium_CrossChannel,
    redeclare package Medium_PipeChannel = Medium_PipeChannel,
    init_pipe=ThermofluidStream.Processes.Internal.InitializationMethodsCondElement.h,
    h_0_pipe=PipeInflowEnthalpy.k,
    T_0_cross=CrossInflowTemperature.k,
    Tin(T(fixed=true))) annotation (Placement(transformation(extent={{-20,-40},{20,0}})));
  ThermofluidStream.Boundaries.Source source_CrossChannel(
    redeclare package Medium = Medium_CrossChannel,
    temperatureFromInput=true,
    T0_par=293.15,
    p0_par=100010) annotation (Placement(transformation(
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
        origin={-30,-60})));
  ThermofluidStream.Boundaries.Sink sink_CrossChannel(redeclare package Medium = Medium_CrossChannel, p0_par=100000) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={82,-20})));
  ThermofluidStream.Boundaries.Sink sink_PipeChannel(redeclare package Medium = Medium_PipeChannel,
    pressureFromInput=true,
    p0_par=100000)                                                                                                 annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={20,90})));
  ThermofluidStream.Processes.FlowResistance flowResistance_Cross(
    redeclare package Medium = Medium_CrossChannel,
    r=0.07,
    l=0.2,
    redeclare function pLoss = ThermofluidStream.Processes.Internal.FlowResistance.laminarTurbulentPressureLossHaaland) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={50,-20})));
  ThermofluidStream.Processes.FlowResistance flowResistance_Pipe(
    redeclare package Medium = Medium_PipeChannel,
    r=0.002,
    l=4,
    redeclare function pLoss = ThermofluidStream.Processes.Internal.FlowResistance.laminarPressureLoss) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=90,
        origin={0,70})));
  inner ThermofluidStream.DropOfCommons dropOfCommons(assertionLevel=AssertionLevel.warning)
                                                                                annotation (Placement(transformation(extent={{60,-100},{80,-80}})));
  ThermofluidStream.FlowControl.MCV mCV_cross(
    redeclare package Medium = Medium_CrossChannel,
    initM_flow=ThermofluidStream.Utilities.Types.InitializationMethods.state,
    m_flow_0=0,
    mode=ThermofluidStream.FlowControl.Internal.Types.MassflowControlValveMode.volume_flow,
    setpointFromInput=true,
    volumeFlow_set_par=0.08333) annotation (Placement(transformation(extent={{-60,-30},{-40,-10}})));
  Modelica.Blocks.Sources.Constant CrossFlowRate(k(
      quantity="VolumeFlowRate",
      unit="m3/s",
      displayUnit="m3/h") = 0.016666666666667,  y(
      quantity="VolumeFlowRate",
      unit="m3/s",
      displayUnit="m3/h")) annotation (Placement(transformation(extent={{-90,30},{-70,50}})));
  Modelica.Blocks.Sources.Constant PipeInflowEnthalpy(k(
      quantity="SpecificEnthalpy",
      unit="J/kg",
      displayUnit="kJ/kg") = 500000, y(
      quantity="SpecificEnthalpy",
      unit="J/kg",
      displayUnit="kJ/kg")) annotation (Placement(transformation(extent={{-92,-94},{-72,-74}})));
  Modelica.Blocks.Sources.Constant PipeInflowPressure(k(
      quantity="AbsolutePressure",
      unit="Pa",
      displayUnit="bar") = 400000,  y(
      quantity="AbsolutePressure",
      unit="Pa",
      displayUnit="bar")) annotation (Placement(transformation(extent={{-92,-60},{-72,-40}})));
  Modelica.Blocks.Sources.Constant PipeMassFlowRate(k(
      quantity="MassFlowRate",
      unit="kg/s",
      displayUnit="kg/h") = 0.00055555555555556,y(
      quantity="MassFlowRate",
      unit="kg/s",
      displayUnit="kg/h")) annotation (Placement(transformation(extent={{-40,34},{-20,54}})));
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
        origin={0,44})));
  Modelica.Blocks.Sources.Constant CrossInflowTemperature(k(
      quantity="Temperature",
      unit="K",
      displayUnit="degC") = 296.35, y(
      quantity="Temperature",
      unit="K",
      displayUnit="degC")) annotation (Placement(transformation(extent={{-70,0},{-90,20}})));
  ThermofluidStream.Sensors.TwoPhaseSensorSelect VaporQuality_Source_Sensor(
    redeclare package Medium = Medium_PipeChannel,
    digits=3,
    quantity=ThermofluidStream.Sensors.Internal.Types.TwoPhaseQuantities.x_kgpkg) annotation (Placement(transformation(extent={{10,-70},{30,-50}})));
  ThermofluidStream.Sensors.TwoPhaseSensorSelect VaporQuality_PipeOut_Sensor(
    redeclare package Medium = Medium_PipeChannel,
    digits=3,
    quantity=ThermofluidStream.Sensors.Internal.Types.TwoPhaseQuantities.x_kgpkg) annotation (Placement(transformation(extent={{-10,0},{-30,20}})));
  Modelica.Blocks.Sources.Constant PipeOutflowPressure(k(
      quantity="AbsolutePressure",
      unit="Pa",
      displayUnit="bar") = 390000, y(
      quantity="AbsolutePressure",
      unit="Pa",
      displayUnit="bar")) annotation (Placement(transformation(extent={{60,80},{40,100}})));
  ThermofluidStream.Sensors.SingleSensorSelect Temperature_PipeIn_Sensor(redeclare package Medium = Medium_PipeChannel, quantity=ThermofluidStream.Sensors.Internal.Types.Quantities.T_C) annotation (Placement(transformation(extent={{10,-90},{30,-70}})));
  ThermofluidStream.Sensors.TwoPhaseSensorSelect TemperatureSaturation_PipeIn_Sensor(redeclare package Medium = Medium_PipeChannel, quantity=ThermofluidStream.Sensors.Internal.Types.TwoPhaseQuantities.T_sat_C) annotation (Placement(transformation(extent={{10,-102},{30,-82}})));
  ThermofluidStream.Sensors.SingleSensorSelect Temperature_PipeOut_Sensor(redeclare package Medium = Medium_PipeChannel, quantity=ThermofluidStream.Sensors.Internal.Types.Quantities.T_C)
                                                                                                                                                                                          annotation (Placement(transformation(extent={{10,0},{30,20}})));
  ThermofluidStream.Sensors.SingleSensorSelect Enthalpy_PipeOut_Sensor(
    redeclare package Medium = Medium_PipeChannel,
    digits=0,
    quantity=ThermofluidStream.Sensors.Internal.Types.Quantities.h_Jpkg) annotation (Placement(transformation(extent={{10,10},{36,30}})));
  ThermofluidStream.Sensors.SingleSensorSelect Temperature_PipeOut_Sensor1(
    redeclare package Medium = Medium_PipeChannel,
    digits=0,
    quantity=ThermofluidStream.Sensors.Internal.Types.Quantities.p_Pa)                                                                                                                    annotation (Placement(transformation(extent={{-8,10},{-40,30}})));
equation
  connect(flowResistance_Pipe.outlet, sink_PipeChannel.inlet) annotation (Line(
      points={{5.55112e-16,80},{4,80},{4,90},{10,90}},
      color={28,108,200},
      thickness=0.5));
  connect(source_PipeChannel.outlet, hXcellPair_Refrig.inlet_Pipe) annotation (Line(
      points={{-20,-60},{0,-60},{0,-40}},
      color={28,108,200},
      thickness=0.5));
  connect(hXcellPair_Refrig.outlet_Cross, flowResistance_Cross.inlet) annotation (Line(
      points={{20,-20},{40,-20}},
      color={28,108,200},
      thickness=0.5));
  connect(flowResistance_Cross.outlet, sink_CrossChannel.inlet) annotation (Line(
      points={{60,-20},{72,-20}},
      color={28,108,200},
      thickness=0.5));
  connect(CrossFlowRate.y, mCV_cross.setpoint_var) annotation (Line(points={{-69,40},{-50,40},{-50,-12}}, color={0,0,127}));
  connect(mCV_cross.outlet, hXcellPair_Refrig.inlet_Cross) annotation (Line(
      points={{-40,-20},{-20,-20}},
      color={28,108,200},
      thickness=0.5));
  connect(source_CrossChannel.outlet, mCV_cross.inlet) annotation (Line(
      points={{-72,-20},{-60,-20}},
      color={28,108,200},
      thickness=0.5));
  connect(source_PipeChannel.p0_var, PipeInflowPressure.y) annotation (Line(points={{-32,-54},{-50,-54},{-50,-50},{-71,-50}}, color={0,0,127}));
  connect(PipeMassFlowRate.y,mCV_pipe. setpoint_var) annotation (Line(points={{-19,44},{-8,44}}, color={0,0,127}));
  connect(hXcellPair_Refrig.outlet_Pipe, mCV_pipe.inlet) annotation (Line(
      points={{0,0},{0,10},{-5.55112e-16,10},{-5.55112e-16,34}},
      color={28,108,200},
      thickness=0.5));
  connect(mCV_pipe.outlet, flowResistance_Pipe.inlet) annotation (Line(
      points={{5.55112e-16,54},{0,54},{0,60},{-5.55112e-16,60}},
      color={28,108,200},
      thickness=0.5));
  connect(source_CrossChannel.T0_var, CrossInflowTemperature.y) annotation (Line(points={{-84,-20},{-100,-20},{-100,10},{-91,10}}, color={0,0,127}));
  connect(PipeInflowEnthalpy.y, source_PipeChannel.h0_var) annotation (Line(points={{-71,-84},{-48,-84},{-48,-60},{-32,-60}}, color={0,0,127}));
  connect(VaporQuality_Source_Sensor.inlet, source_PipeChannel.outlet) annotation (Line(
      points={{10,-60},{-20,-60}},
      color={28,108,200},
      thickness=0.5));
  connect(hXcellPair_Refrig.outlet_Pipe, VaporQuality_PipeOut_Sensor.inlet) annotation (Line(
      points={{0,0},{0,10},{-10,10}},
      color={28,108,200},
      thickness=0.5));
  connect(sink_PipeChannel.p0_var, PipeOutflowPressure.y) annotation (Line(points={{22,90},{39,90}}, color={0,0,127}));
  connect(TemperatureSaturation_PipeIn_Sensor.inlet, source_PipeChannel.outlet) annotation (Line(
      points={{10,-92},{0,-92},{0,-60},{-20,-60}},
      color={28,108,200},
      thickness=0.5));
  connect(Temperature_PipeIn_Sensor.inlet, source_PipeChannel.outlet) annotation (Line(
      points={{10,-80},{0,-80},{0,-60},{-20,-60}},
      color={28,108,200},
      thickness=0.5));
  connect(hXcellPair_Refrig.outlet_Pipe, Temperature_PipeOut_Sensor.inlet) annotation (Line(
      points={{0,0},{0,10},{10,10}},
      color={28,108,200},
      thickness=0.5));
  connect(Enthalpy_PipeOut_Sensor.inlet, hXcellPair_Refrig.outlet_Pipe) annotation (Line(
      points={{10,20},{0,20},{0,0}},
      color={28,108,200},
      thickness=0.5));
  connect(hXcellPair_Refrig.outlet_Pipe, Temperature_PipeOut_Sensor1.inlet) annotation (Line(
      points={{0,0},{0,20},{-8,20}},
      color={28,108,200},
      thickness=0.5));
  annotation (
    Icon(coordinateSystem(preserveAspectRatio=false)),
    Diagram(coordinateSystem(preserveAspectRatio=false), graphics={
        Text(
          extent={{26,72},{66,50}},
          textColor={28,108,200},
          textString="Pipe Flow"),
        Text(
          extent={{40,-40},{80,-62}},
          textColor={28,108,200},
          textString="Cross Flow"),
        Text(
          extent={{-66,-90},{-34,-96}},
          textColor={28,108,200},
          textString="2-phase region:
10 bar: 250 .. 420 kJ/kg"),
        Text(
          extent={{-42,70},{-14,58}},
          textColor={28,108,200},
          textString="GUNT: 10.5 kg/s")}),
    experiment(
      StopTime=100,
      Interval=0.01,
      __Dymola_Algorithm="Esdirk23a"),
    __Dymola_Commands(file="Scripts/Plots_HXcellPair_Refrig.mos"));
end Test_HXcellPair_Condensor;
