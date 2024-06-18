within ThermofluidStream.HeatExchangersPhysical.HXutilities.Controller;
model Test_Controller
    extends Modelica.Icons.Example;
  PID_varLimits PID(
    controllerType=Controller.ControllerType.PID,
    useMaxPort=true,
    k=5,
    Ti=10) annotation (Placement(transformation(extent={{-20,-20},{20,20}})));
  Modelica.Blocks.Sources.Step           step(
    height=1,
    offset=0,
    startTime=100)
    annotation (Placement(transformation(extent={{-68,-10},{-48,10}})));
  Modelica.Blocks.Sources.RealExpression MeasuredValue(y=0)
    annotation (Placement(transformation(extent={{-68,-60},{-48,-40}})));
  Modelica.Blocks.Interaction.Show.RealValue Actuating_Signal(significantDigits=5)
    annotation (Placement(transformation(extent={{50,-12},{80,12}})));
  Modelica.Blocks.Sources.RealExpression UpperLimit(y=10)
    annotation (Placement(transformation(extent={{70,6},{50,26}})));
equation
  connect(PID.y, Actuating_Signal.numberPort) annotation (Line(points={{22,0},{47.75,0}}, color={0,0,127}));
  connect(UpperLimit.y, PID.y_max) annotation (Line(points={{49,16},{22,16}}, color={0,0,127}));
  connect(step.y, PID.u_s) annotation (Line(points={{-47,0},{-22,0}}, color={0,0,127}));
  connect(MeasuredValue.y, PID.u_m) annotation (Line(points={{-47,-50},{0,-50},{0,-22}}, color={0,0,127}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)),
    experiment(
      StopTime=200,
      Interval=0.5,
      __Dymola_Algorithm="Dassl"));
end Test_Controller;
