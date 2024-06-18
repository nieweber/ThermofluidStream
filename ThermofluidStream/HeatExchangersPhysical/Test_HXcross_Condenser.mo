within ThermofluidStream.HeatExchangersPhysical;
model Test_HXcross_Condenser "Test of cross heat exchanger moist air to refrigerant"
    extends Modelica.Icons.Example;

  replaceable package Medium_CrossChannel =
      ThermofluidStream.Media.myMedia.Air.MoistAir                        constrainedby ThermofluidStream.Media.myMedia.Interfaces.PartialMedium
                                                        annotation(choicesAllMatching = true);

//  package Medium_PipeChannel = TILMediaWrapper.VLEFluidWrapper.R513A;
  replaceable package Medium_PipeChannel = ThermofluidStream.Media.myMedia.R134a.R134a_ph
                                                                          constrainedby ThermofluidStream.Media.myMedia.Interfaces.PartialTwoPhaseMedium
                                                        annotation(choicesAllMatching = true);

  HX_1PhaseCross_2PhasePipe hXcross_Refrig(
    redeclare package Medium_CrossChannel = Medium_CrossChannel,
    redeclare package Medium_PipeChannel = Medium_PipeChannel,
    init_pipe=ThermofluidStream.Processes.Internal.InitializationMethodsCondElement.h,
    h_0_pipe=PipeInflowEnthalpy.k,
    T_0_cross=CrossInflowTemperature.k) annotation (Placement(transformation(extent={{-20,-40},{20,0}})));
  ThermofluidStream.Boundaries.Source source_CrossChannel(
    redeclare package Medium = Medium_CrossChannel,
    temperatureFromInput=true,
    T0_par=293.15,
    p0_par=102000) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={-100,-20})));
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
        origin={-40,-70})));
  ThermofluidStream.Boundaries.Sink sink_CrossChannel(redeclare package Medium = Medium_CrossChannel, p0_par=101300) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={110,-20})));
  ThermofluidStream.Boundaries.Sink sink_PipeChannel(redeclare package Medium = Medium_PipeChannel, p0_par=390000) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=90,
        origin={0,100})));
  ThermofluidStream.Processes.FlowResistance flowResistance_Cross(
    redeclare package Medium = Medium_CrossChannel,
    r=0.07,
    l=0.2,
    redeclare function pLoss = ThermofluidStream.Processes.Internal.FlowResistance.laminarTurbulentPressureLossHaaland) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={80,-20})));
  ThermofluidStream.Processes.FlowResistance flowResistance_Pipe(
    redeclare package Medium = Medium_PipeChannel,
    initM_flow=ThermofluidStream.Utilities.Types.InitializationMethods.state,
    r=0.002,
    l=4,
    redeclare function pLoss = ThermofluidStream.Processes.Internal.FlowResistance.laminarPressureLoss) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=90,
        origin={0,74})));
  inner ThermofluidStream.DropOfCommons dropOfCommons(assertionLevel=AssertionLevel.warning)
                                                                                annotation (Placement(transformation(extent={{60,-100},{80,-80}})));
  ThermofluidStream.FlowControl.MCV mCV_cross(
    redeclare package Medium = Medium_CrossChannel,
    mode=ThermofluidStream.FlowControl.Internal.Types.MassflowControlValveMode.volume_flow,
    setpointFromInput=true,
    volumeFlow_set_par=0.08333) annotation (Placement(transformation(extent={{-80,-30},{-60,-10}})));
  ThermofluidStream.FlowControl.MCV mCV_pipe(
    redeclare package Medium = Medium_PipeChannel,
    mode=ThermofluidStream.FlowControl.Internal.Types.MassflowControlValveMode.mass_flow,
    setpointFromInput=true,
    volumeFlow_set_par(displayUnit="l/h") = 1.9444444444444e-06,
    dp_int(start=-1000))                                         annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=90,
        origin={0,50})));
  Modelica.Blocks.Sources.Constant PipeFlowRate(k(
      quantity="VolumeFlowRate",
      unit="m3/s",
      displayUnit="l/h") = 1.9444444444444e-06,
                            y(
      quantity="VolumeFlowRate",
      unit="m3/s",
      displayUnit="l/h")) annotation (Placement(transformation(extent={{-80,70},{-60,90}})));
  Modelica.Blocks.Sources.Constant CrossFlowRate(k(
      quantity="VolumeFlowRate",
      unit="m3/s",
      displayUnit="m3/h") = 0.13333333333333,
                                 y(
      quantity="VolumeFlowRate",
      unit="m3/s",
      displayUnit="m3/h")) annotation (Placement(transformation(extent={{-100,40},{-80,60}})));
  Modelica.Blocks.Sources.Constant PipeInflowPressure(k(
      quantity="AbsolutePressure",
      unit="Pa",
      displayUnit="bar") = 800000,
                                 y(
      quantity="AbsolutePressure",
      unit="Pa",
      displayUnit="bar")) annotation (Placement(transformation(extent={{-100,-60},{-80,-40}})));
  Modelica.Blocks.Sources.Constant PipeMassFlowRate(k(
      quantity="MassFlowRate",
      unit="kg/s",
      displayUnit="kg/h") = 0.0027777777777778,  y(
      quantity="MassFlowRate",
      unit="kg/s",
      displayUnit="kg/h")) annotation (Placement(transformation(extent={{-40,40},{-20,60}})));
  Modelica.Blocks.Interaction.Show.RealValue MassFlowValue(use_numberPort=false, number(
      quantity="MassFlowRate",
      unit="kg/s",
      displayUnit="kg/h")=mCV_pipe.m_flow,
    significantDigits=3)                   annotation (Placement(transformation(extent={{20,54},{40,70}})));
  Modelica.Blocks.Interaction.Show.RealValue VolumeFlowValue(use_numberPort=false, number(
      quantity="VolumeFlowRate",
      unit="m3/s",
      displayUnit="l/h")=mCV_pipe.V_flow,
    significantDigits=3)                  annotation (Placement(transformation(extent={{20,42},{40,58}})));
  Modelica.Blocks.Sources.Constant CrossInflowTemperature(k(
      quantity="Temperature",
      unit="K",
      displayUnit="degC") = 296.35, y(
      quantity="Temperature",
      unit="K",
      displayUnit="degC")) annotation (Placement(transformation(extent={{-140,-30},{-120,-10}})));
  Modelica.Blocks.Sources.Constant PipeInflowEnthalpy(k(
      quantity="SpecificEnthalpy",
      unit="J/kg",
      displayUnit="kJ/kg") = 450000, y(
      quantity="SpecificEnthalpy",
      unit="J/kg",
      displayUnit="kJ/kg")) annotation (Placement(transformation(extent={{-100,-94},{-80,-74}})));
  ThermofluidStream.Sensors.SingleSensorSelect Temperature_CrossOut_Sensor(redeclare package Medium = Medium_CrossChannel, quantity=ThermofluidStream.Sensors.Internal.Types.Quantities.T_C) annotation (Placement(transformation(extent={{40,-40},{60,-20}})));
  ThermofluidStream.Sensors.SingleFlowSensor EnthalpyFlow_CrossIn_Sensor(
    redeclare package Medium = Medium_CrossChannel,
    digits=0,
    quantity=ThermofluidStream.Sensors.Internal.Types.MassFlowQuantities.H_flow_Jps) annotation (Placement(transformation(extent={{-50,-24},{-30,-4}})));
  ThermofluidStream.Sensors.SingleFlowSensor EnthalpyFlow_CrossOut_Sensor(
    redeclare package Medium = Medium_CrossChannel,
    digits=0,
    quantity=ThermofluidStream.Sensors.Internal.Types.MassFlowQuantities.H_flow_Jps) annotation (Placement(transformation(extent={{40,-24},{60,-4}})));
  ThermofluidStream.Sensors.SingleFlowSensor EnthalpyFlow_PipeOut_Sensor(
    redeclare package Medium = Medium_PipeChannel,
    digits=0,
    quantity=ThermofluidStream.Sensors.Internal.Types.MassFlowQuantities.H_flow_Jps) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=90,
        origin={-6,24})));
  ThermofluidStream.Sensors.SingleFlowSensor EnthalpyFlow_PipeIn_Sensor(
    redeclare package Medium = Medium_PipeChannel,
    digits=0,
    quantity=ThermofluidStream.Sensors.Internal.Types.MassFlowQuantities.H_flow_Jps) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=90,
        origin={-6,-56})));
  ThermofluidStream.Sensors.SingleSensorSelect Temperature_PipeOut_Sensor(redeclare package Medium = Medium_PipeChannel, quantity=ThermofluidStream.Sensors.Internal.Types.Quantities.T_C) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={20,8})));
  ThermofluidStream.Sensors.SingleSensorSelect Temperature_PipeIn_Sensor(redeclare package Medium = Medium_PipeChannel, quantity=ThermofluidStream.Sensors.Internal.Types.Quantities.T_C) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={20,-70})));
  ThermofluidStream.Sensors.TwoPhaseSensorSelect TemperatureSaturation_PipeIn_Sensor(redeclare package Medium = Medium_PipeChannel, quantity=ThermofluidStream.Sensors.Internal.Types.TwoPhaseQuantities.T_sat_C) annotation (Placement(transformation(extent={{10,-90},{30,-70}})));
  ThermofluidStream.Sensors.TwoPhaseSensorSelect VaporQuality_Source_Sensor(
    redeclare package Medium = Medium_PipeChannel,
    digits=3,
    quantity=ThermofluidStream.Sensors.Internal.Types.TwoPhaseQuantities.x_kgpkg) annotation (Placement(transformation(extent={{-24,-90},{-4,-70}})));
  ThermofluidStream.Sensors.TwoPhaseSensorSelect VaporQuality_PipeOut_Sensor(
    redeclare package Medium = Medium_PipeChannel,
    digits=3,
    quantity=ThermofluidStream.Sensors.Internal.Types.TwoPhaseQuantities.x_kgpkg) annotation (Placement(transformation(extent={{-10,-2},{-30,18}})));
equation
  connect(flowResistance_Pipe.outlet, sink_PipeChannel.inlet) annotation (Line(
      points={{0,84},{0,90}},
      color={28,108,200},
      thickness=0.5));
  connect(flowResistance_Cross.outlet, sink_CrossChannel.inlet) annotation (Line(
      points={{90,-20},{100,-20}},
      color={28,108,200},
      thickness=0.5));
  connect(source_CrossChannel.outlet, mCV_cross.inlet) annotation (Line(
      points={{-90,-20},{-80,-20}},
      color={28,108,200},
      thickness=0.5));
  connect(flowResistance_Pipe.inlet, mCV_pipe.outlet) annotation (Line(
      points={{0,64},{0,60}},
      color={28,108,200},
      thickness=0.5));
  connect(CrossFlowRate.y, mCV_cross.setpoint_var) annotation (Line(points={{-79,50},{-70,50},{-70,-12}}, color={0,0,127}));
  connect(PipeInflowPressure.y, source_PipeChannel.p0_var) annotation (Line(points={{-79,-50},{-60,-50},{-60,-64},{-42,-64}}, color={0,0,127}));
  connect(PipeMassFlowRate.y, mCV_pipe.setpoint_var) annotation (Line(points={{-19,50},{-8,50}}, color={0,0,127}));
  connect(CrossInflowTemperature.y, source_CrossChannel.T0_var) annotation (Line(points={{-119,-20},{-102,-20}},                 color={0,0,127}));
  connect(PipeInflowEnthalpy.y, source_PipeChannel.h0_var) annotation (Line(points={{-79,-84},{-60,-84},{-60,-70},{-42,-70}}, color={0,0,127}));
  connect(Temperature_CrossOut_Sensor.inlet, hXcross_Refrig.outlet_Cross) annotation (Line(
      points={{40,-30},{30,-30},{30,-20},{20.8,-20}},
      color={28,108,200},
      thickness=0.5));
  connect(EnthalpyFlow_CrossIn_Sensor.outlet, hXcross_Refrig.inlet_Cross) annotation (Line(
      points={{-30,-14},{-26,-14},{-26,-20},{-20,-20}},
      color={28,108,200},
      thickness=0.5));
  connect(EnthalpyFlow_CrossIn_Sensor.inlet, mCV_cross.outlet) annotation (Line(
      points={{-50,-14},{-56,-14},{-56,-20},{-60,-20}},
      color={28,108,200},
      thickness=0.5));
  connect(EnthalpyFlow_CrossOut_Sensor.outlet, flowResistance_Cross.inlet) annotation (Line(
      points={{60,-14},{66,-14},{66,-20},{70,-20}},
      color={28,108,200},
      thickness=0.5));
  connect(EnthalpyFlow_CrossOut_Sensor.inlet, hXcross_Refrig.outlet_Cross) annotation (Line(
      points={{40,-14},{30,-14},{30,-20},{20.8,-20}},
      color={28,108,200},
      thickness=0.5));
  connect(EnthalpyFlow_PipeOut_Sensor.inlet, hXcross_Refrig.outlet_Pipe) annotation (Line(
      points={{-6,14},{-6,3.2},{0,3.2},{0,0.4}},
      color={28,108,200},
      thickness=0.5));
  connect(EnthalpyFlow_PipeOut_Sensor.outlet, mCV_pipe.inlet) annotation (Line(
      points={{-6,34},{-6,36},{-5.55112e-16,36},{-5.55112e-16,40}},
      color={28,108,200},
      thickness=0.5));
  connect(EnthalpyFlow_PipeIn_Sensor.inlet, source_PipeChannel.outlet) annotation (Line(
      points={{-6,-66},{0,-66},{0,-70},{-30,-70}},
      color={28,108,200},
      thickness=0.5));
  connect(EnthalpyFlow_PipeIn_Sensor.outlet, hXcross_Refrig.inlet_Pipe) annotation (Line(
      points={{-6,-46},{-6,-43.2},{0,-43.2},{0,-40.4}},
      color={28,108,200},
      thickness=0.5));
  connect(hXcross_Refrig.outlet_Pipe, Temperature_PipeOut_Sensor.inlet) annotation (Line(
      points={{0,0.4},{0,8},{10,8}},
      color={28,108,200},
      thickness=0.5));
  connect(Temperature_PipeIn_Sensor.inlet, source_PipeChannel.outlet) annotation (Line(
      points={{10,-70},{-30,-70}},
      color={28,108,200},
      thickness=0.5));
  connect(TemperatureSaturation_PipeIn_Sensor.inlet, source_PipeChannel.outlet) annotation (Line(
      points={{10,-80},{0,-80},{0,-70},{-30,-70}},
      color={28,108,200},
      thickness=0.5));
  connect(VaporQuality_Source_Sensor.inlet, source_PipeChannel.outlet) annotation (Line(
      points={{-24,-80},{-26,-80},{-26,-70},{-30,-70}},
      color={28,108,200},
      thickness=0.5));
  connect(VaporQuality_PipeOut_Sensor.inlet, hXcross_Refrig.outlet_Pipe) annotation (Line(
      points={{-10,8},{0,8},{0,0.4}},
      color={28,108,200},
      thickness=0.5));
  annotation (
    Icon(coordinateSystem(preserveAspectRatio=false, extent={{-160,-100},{120,100}})),
    Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-160,-100},{120,100}}), graphics={
        Text(
          extent={{64,98},{104,76}},
          textColor={0,0,0},
          textString="Pipe Flow"),
        Text(
          extent={{74,-32},{114,-54}},
          textColor={0,0,0},
          textString="Cross Flow"),
        Text(
          extent={{-74,-90},{-42,-96}},
          textColor={28,108,200},
          textString="2-phase region:
250 .. 420 kJ/kg"),
        Text(
          extent={{-42,76},{-14,64}},
          textColor={28,108,200},
          textString="GUNT: 10.5 kg/s"),
        Text(
          extent={{-74,-104},{-42,-110}},
          textColor={28,108,200},
          textString="Condenser Input:
430.8 kJ/kg
10 bar"),
        Text(
          extent={{26,114},{58,108}},
          textColor={28,108,200},
          textString="Condenser Input:
9.995 bar"),
        Text(
          extent={{-108,2},{-84,-12}},
          textColor={0,0,0},
          textString="Moist Air"),
        Text(
          extent={{-50,-48},{-26,-62}},
          textColor={0,0,0},
          textString="Refrigerant"),
        Text(
          extent={{98,0},{116,-12}},
          textColor={0,0,0},
          textString="Moist Air"),
        Text(
          extent={{-36,104},{-12,90}},
          textColor={0,0,0},
          textString="Refrigerant")}),
    experiment(StopTime=1000, __Dymola_Algorithm="Dassl"),
    __Dymola_Commands(file="Scripts/Plots_HXcross_Refrig.mos"));
end Test_HXcross_Condenser;
