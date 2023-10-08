within TFS_Tutorial.Submodels;
model VCS

  package Coolant =
      ThermofluidStream.Media.myMedia.Incompressible.Examples.Glycol47;
  package Air = ThermofluidStream.Media.myMedia.Air.MoistAir;
  package Refrigerant = ThermofluidStream.Media.XRGMedia.R134a_ph;

  parameter Modelica.Units.SI.Time T_PI = 1;
  parameter Real k_PI(unit="1") = 1;

  ThermofluidStream.HeatExchangers.DiscretizedCounterFlowHEX condenser(
    redeclare package MediumA = Air,
    redeclare package MediumB = Refrigerant,
    redeclare model ConductionElementB = ThermofluidStream.HeatExchangers.Internal.ConductionElementHEX_twoPhase (U_tp_nom=3000, m_flow_nom=0.1),
    initializeMassFlow=false,
    nCells=5,
    A=30,
    k_wall=100)  annotation (Placement(transformation(
        extent={{10,10},{-10,-10}},
        rotation=0,
        origin={0,112})));
  ThermofluidStream.Processes.Compressor compressor(
    redeclare package Medium = Refrigerant,
    initM_flow=ThermofluidStream.Utilities.Types.InitializationMethods.none,
    omega_from_input=true,
    redeclare function dp_tau_compressor = ThermofluidStream.Processes.Internal.TurboComponent.dp_tau_const_isentrop (
        eta=0.7,
        kappaFromMedia=false,
        kappa_fixed=1.13))
    annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=90,
        origin={68,0})));
  ThermofluidStream.Examples.Utilities.Receiver receiver(
    redeclare package Medium = Refrigerant,
    p_start=1000000,
    V_par=0.001,
    init_method=ThermofluidStream.Boundaries.Internal.InitializationMethodsPhaseSeperator.l)
    annotation (Placement(transformation(
        extent={{-10,10},{10,-10}},
        rotation=180,
        origin={-52,104})));
  ThermofluidStream.FlowControl.BasicControlValve basicControlValve(
    redeclare package Medium = Refrigerant,
    initM_flow=ThermofluidStream.Utilities.Types.InitializationMethods.state,
    m_flow_0=0.3,
    invertInput=true,
    redeclare function valveCharacteristics = ThermofluidStream.FlowControl.Internal.ControlValve.linearCharacteristics,
    flowCoefficient=ThermofluidStream.FlowControl.Internal.Types.FlowCoefficientTypesBasic.m_flow_set,
    Kvs(displayUnit="m3/s"),
    m_flow_ref_set=0.2)  annotation (Placement(transformation(
        extent={{10,10},{-10,-10}},
        rotation=90,
        origin={-104,0})));
  ThermofluidStream.HeatExchangers.DiscretizedCounterFlowHEX evaporator(
    redeclare package MediumA = Coolant,
    redeclare package MediumB = Refrigerant,
    redeclare model ConductionElementB = ThermofluidStream.HeatExchangers.Internal.ConductionElementHEX_twoPhase (U_tp_nom=2000, m_flow_nom=0.1),
    initializeMassFlow=false,
    nCells=5,
    A=1,
    V_Hex(displayUnit="l"),
    k_wall=300) annotation (Placement(transformation(
        extent={{10,10},{-10,-10}},
        rotation=180,
        origin={-2,-112})));
  ThermofluidStream.Sensors.SingleSensorSelect
                               sensorVaporQuality11(
    redeclare package Medium = Refrigerant,
    outputValue=true,
    quantity=ThermofluidStream.Sensors.Internal.Types.Quantities.p_bar,
    filter_output=true) annotation (Placement(transformation(extent={{10,10},{-10,-10}},
        rotation=180,
        origin={90,42})));
  ThermofluidStream.Sensors.TwoPhaseSensorSelect
                               sensorVaporQuality8(
    redeclare package Medium = Refrigerant,
    outputValue=true,
    quantity=ThermofluidStream.Sensors.Internal.Types.TwoPhaseQuantities.T_oversat_K,
    filter_output=true) annotation (Placement(transformation(extent={{10,10},{-10,-10}},
        rotation=270,
        origin={14,-70})));
  Modelica.Blocks.Continuous.LimPID PI(
    controllerType=Modelica.Blocks.Types.SimpleController.PI,
    k=k_PI,
    Ti=T_PI,
    yMax=1,
    yMin=0,
    initType=Modelica.Blocks.Types.Init.InitialOutput,
    y_start=0.15)
               annotation (Placement(transformation(extent={{-16,-6},{-28,6}})));
  Modelica.Blocks.Sources.Constant const(k=5) annotation (Placement(transformation(extent={{6,-6},{-6,6}})));
  ThermofluidStream.Interfaces.Inlet inletAirCondenser(redeclare package Medium = Air) annotation (Placement(transformation(extent={{-168,110},{-148,130}}), iconTransformation(extent={{-168,110},{-148,130}})));
  ThermofluidStream.Interfaces.Inlet inletAirEvaporator(redeclare package Medium = Coolant)
                                                                                        annotation (Placement(transformation(extent={{168,-130},{148,-110}}), iconTransformation(extent={{168,-130},{148,-110}})));
  ThermofluidStream.Interfaces.Outlet outletAirEvaporator(redeclare package Medium = Coolant)
                                                                                          annotation (Placement(transformation(extent={{-152,-130},{-172,-110}}), iconTransformation(extent={{-152,-130},{-172,-110}})));
  ThermofluidStream.Interfaces.Outlet outletAirCondenser(redeclare package Medium = Air) annotation (Placement(transformation(extent={{150,110},{170,130}}), iconTransformation(extent={{150,110},{170,130}})));
  Controllers.LimPID             PID(
    k=-100,
    Ti=T_PI*0.1,
    Td=1000*T_PI,
    yMax=10000,
    yMin=100,
    initType=Modelica.Blocks.Types.Init.InitialOutput,
    y_start=1500)                                      annotation (Placement(
        transformation(
        extent={{-6,6},{6,-6}},
        rotation=180,
        origin={108,0})));
  Modelica.Blocks.Interfaces.RealInput u_m annotation (Placement(transformation(extent={{168,-44},{142,-18}}), iconTransformation(extent={{168,-44},{142,-18}})));
  Modelica.Blocks.Interfaces.RealInput u_set annotation (Placement(transformation(extent={{168,16},{142,42}}), iconTransformation(extent={{168,16},{142,42}})));
  ThermofluidStream.Sensors.MultiSensor_Tpm multiSensor_Tpm6(
    redeclare package Medium = Refrigerant,
    temperatureUnit="degC",
    pressureUnit="bar")
    annotation (Placement(transformation(extent={{-10,-10},{10,10}},
        rotation=180,
        origin={46,94})));
  ThermofluidStream.Sensors.SingleSensorSelect
                                            singleSensorSelect(redeclare package Medium = Refrigerant, quantity=ThermofluidStream.Sensors.Internal.Types.Quantities.T_C)
    annotation (Placement(transformation(extent={{-10,10},{10,-10}},
        rotation=180,
        origin={-36,-92})));
  ThermofluidStream.Sensors.SingleSensorSelect
                                            singleSensorSelect1(redeclare package Medium = Refrigerant, quantity=ThermofluidStream.Sensors.Internal.Types.Quantities.T_C)
    annotation (Placement(transformation(extent={{10,10},{-10,-10}},
        rotation=180,
        origin={40,-88})));
  ThermofluidStream.Sensors.MultiSensor_Tpm multiSensor_Tpm4(
    redeclare package Medium = Coolant,
    temperatureUnit="degC",
    pressureUnit="bar") annotation (Placement(transformation(
        extent={{10,10},{-10,-10}},
        rotation=0,
        origin={92,-130})));
  ThermofluidStream.Sensors.MultiSensor_Tpm multiSensor_Tpm1(
    redeclare package Medium = Coolant,
    temperatureUnit="degC",
    pressureUnit="bar") annotation (Placement(transformation(
        extent={{10,10},{-10,-10}},
        rotation=0,
        origin={-112,-130})));
  ThermofluidStream.Sensors.MultiSensor_Tpm multiSensor_Tpm2(
    redeclare package Medium = Air,
    temperatureUnit="degC",
    pressureUnit="bar") annotation (Placement(transformation(
        extent={{10,-10},{-10,10}},
        rotation=180,
        origin={118,110})));
  ThermofluidStream.Sensors.MultiSensor_Tpm multiSensor_Tpm3(
    redeclare package Medium = Air,
    temperatureUnit="degC",
    pressureUnit="bar") annotation (Placement(transformation(
        extent={{10,-10},{-10,10}},
        rotation=180,
        origin={-120,110})));
equation
  connect(PI.u_s,const. y) annotation (Line(points={{-14.8,0},{-6.6,0}},                      color={0,0,127}));
  connect(condenser.outletB, receiver.inlet) annotation (Line(
      points={{-10,104},{-42,104}},
      color={28,108,200},
      thickness=0.5));
  connect(receiver.outlet, basicControlValve.inlet) annotation (Line(
      points={{-62,104},{-104,104},{-104,10}},
      color={28,108,200},
      thickness=0.5));
  connect(basicControlValve.outlet, evaporator.inletB) annotation (Line(
      points={{-104,-10},{-104,-104},{-12,-104}},
      color={28,108,200},
      thickness=0.5));
  connect(sensorVaporQuality11.inlet, compressor.outlet) annotation (Line(
      points={{80,42},{68,42},{68,10}},
      color={28,108,200},
      thickness=0.5));
  connect(sensorVaporQuality8.inlet, evaporator.outletB) annotation (Line(
      points={{14,-80},{14,-104},{8,-104}},
      color={28,108,200},
      thickness=0.5));
  connect(PI.u_m, sensorVaporQuality8.value_out) annotation (Line(points={{-22,-7.2},{-22,-44},{14,-44},{14,-60}}, color={0,0,127}));
  connect(compressor.outlet, multiSensor_Tpm6.inlet) annotation (Line(
      points={{68,10},{68,104},{56,104}},
      color={28,108,200},
      thickness=0.5));
  connect(multiSensor_Tpm6.outlet, condenser.inletB) annotation (Line(
      points={{36,104},{10,104}},
      color={28,108,200},
      thickness=0.5));
  connect(singleSensorSelect.inlet, evaporator.inletB) annotation (Line(
      points={{-26,-92},{-20,-92},{-20,-104},{-12,-104}},
      color={28,108,200},
      thickness=0.5));
  connect(singleSensorSelect1.inlet, evaporator.outletB) annotation (Line(
      points={{30,-88},{16,-88},{16,-104},{8,-104}},
      color={28,108,200},
      thickness=0.5));
  connect(evaporator.outletB, compressor.inlet) annotation (Line(
      points={{8,-104},{68,-104},{68,-10}},
      color={28,108,200},
      thickness=0.5));
  connect(inletAirEvaporator, multiSensor_Tpm4.inlet) annotation (Line(
      points={{158,-120},{102,-120}},
      color={28,108,200},
      thickness=0.5));
  connect(multiSensor_Tpm4.outlet, evaporator.inletA) annotation (Line(
      points={{82,-120},{8,-120}},
      color={28,108,200},
      thickness=0.5));
  connect(evaporator.outletA, multiSensor_Tpm1.inlet) annotation (Line(
      points={{-12,-120},{-102,-120}},
      color={28,108,200},
      thickness=0.5));
  connect(multiSensor_Tpm1.outlet, outletAirEvaporator) annotation (Line(
      points={{-122,-120},{-162,-120}},
      color={28,108,200},
      thickness=0.5));
  connect(condenser.outletA, multiSensor_Tpm2.inlet) annotation (Line(
      points={{10,120},{108,120}},
      color={28,108,200},
      thickness=0.5));
  connect(multiSensor_Tpm2.outlet, outletAirCondenser) annotation (Line(
      points={{128,120},{160,120}},
      color={28,108,200},
      thickness=0.5));
  connect(inletAirCondenser, multiSensor_Tpm3.inlet) annotation (Line(
      points={{-158,120},{-130,120}},
      color={28,108,200},
      thickness=0.5));
  connect(multiSensor_Tpm3.outlet, condenser.inletA) annotation (Line(
      points={{-110,120},{-10,120}},
      color={28,108,200},
      thickness=0.5));
  connect(u_m, PID.u_m) annotation (Line(points={{155,-31},{108,-31},{108,-7.2}},  color={0,0,127}));
  connect(u_set, PID.u_s) annotation (Line(points={{155,29},{136,29},{136,-1.11022e-15},{115.2,-1.11022e-15}},
                                                                                             color={0,0,127}));
  connect(PI.y, basicControlValve.u_in) annotation (Line(points={{-28.6,0},{-96,0}}, color={0,0,127}));
  connect(PID.y, compressor.omega_input) annotation (Line(points={{101.4,4.44089e-16},{90,4.44089e-16},{90,-6.66134e-16},{78,-6.66134e-16}},
                                                                                         color={0,0,127}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-160,-160},{160,160}}), graphics={
        Rectangle(
          extent={{-160,160},{160,-160}},
          lineColor={28,108,200},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-40,140},{40,100}},
          lineColor={0,0,0},
          lineThickness=0.5),
        Rectangle(
          extent={{-40,-100},{40,-140}},
          lineColor={0,0,0},
          lineThickness=0.5),
        Ellipse(
          extent={{76,24},{124,-24}},
          lineColor={0,0,0},
          lineThickness=0.5,
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid),
        Line(
          points={{-52,-7},{-14,1}},
          color={0,0,0},
          thickness=0.5,
          origin={111,36},
          rotation=90),
        Polygon(
          points={{-100,0},{-112,-24},{-88,-24},{-100,0}},
          lineColor={0,0,0},
          lineThickness=0.5),
        Polygon(
          points={{-100,0},{-88,24},{-112,24},{-100,0}},
          lineColor={0,0,0},
          lineThickness=0.5),
        Line(
          points={{-100,-24},{-100,-120}},
          color={0,0,0},
          thickness=0.5),
        Line(
          points={{-40,-120},{-100,-120}},
          color={0,0,0},
          thickness=0.5),
        Line(
          points={{100,-120},{40,-120}},
          color={0,0,0},
          thickness=0.5),
        Line(
          points={{100,-24},{100,-120}},
          color={0,0,0},
          thickness=0.5),
        Line(
          points={{100,120},{100,24}},
          color={0,0,0},
          thickness=0.5),
        Line(
          points={{100,120},{40,120}},
          color={0,0,0},
          thickness=0.5),
        Line(
          points={{-40,120},{-100,120}},
          color={0,0,0},
          thickness=0.5),
        Line(
          points={{-100,120},{-100,24}},
          color={0,0,0},
          thickness=0.5),
        Line(
          points={{-52,9},{-14,1}},
          color={0,0,0},
          thickness=0.5,
          origin={91,36},
          rotation=90),
        Text(
          extent={{-34,134},{32,108}},
          textColor={0,0,0},
          textString="COND"),
        Text(
          extent={{-32,-106},{34,-132}},
          textColor={0,0,0},
          textString="EVAP")}),                                  Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-160,-160},{160,160}})));
end VCS;
