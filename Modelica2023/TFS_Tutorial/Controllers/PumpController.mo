within TFS_Tutorial.Controllers;
model PumpController

  parameter Real k_PI = -20 "Gain of PI-Controller";
  parameter Modelica.Units.SI.Time Ti_PI = 5 "Time constant of PI-Controller";

  PIController controller(k=k_PI, Ti=Ti_PI)
                                           annotation (Placement(transformation(extent={{72,-16},{38,16}})));
  Modelica.Blocks.Math.Abs abs1 annotation (Placement(transformation(extent={{-62,-10},{-82,10}})));
  Modelica.Blocks.Logical.Switch switch1 annotation (Placement(transformation(extent={{-28,-10},{-48,10}})));
  Modelica.Blocks.Sources.RealExpression realExpression9(y=5) annotation (Placement(transformation(
        extent={{10,-10},{-10,10}},
        rotation=0,
        origin={0,-20})));
  Modelica.Blocks.Continuous.FirstOrder firstOrder1(T=0.01, initType=Modelica.Blocks.Types.Init.InitialState)
                                                            annotation (Placement(transformation(
        extent={{10,-10},{-10,10}},
        rotation=0,
        origin={-102,0})));
  Modelica.Blocks.Interfaces.RealInput u_m(unit = "K") annotation (Placement(transformation(
        rotation=270,
        extent={{14,-14},{-14,14}},
        origin={-1.77636e-15,-144}),
                          iconTransformation(
        extent={{24,-24},{-24,24}},
        rotation=270,
        origin={0,-136})));
  Modelica.Blocks.Interfaces.RealInput u_s(unit="K") annotation (Placement(transformation(
        rotation=0,
        extent={{24,-24},{-24,24}},
        origin={144,0}), iconTransformation(
        extent={{24,-24},{-24,24}},
        rotation=0,
        origin={136,0})));
  Modelica.Blocks.Interfaces.RealOutput y(start=firstOrder1.y_start) annotation (Placement(transformation(rotation=0, extent={{-132,-18},{-168,18}}), iconTransformation(extent={{-132,-18},{-168,18}})));
  Modelica.Blocks.Interfaces.BooleanInput u_switch annotation (Placement(transformation(
        rotation=90,
        extent={{-24,-24},{24,24}},
        origin={-100,-136}), iconTransformation(
        extent={{-24,-24},{24,24}},
        rotation=90,
        origin={-100,-136})));
equation
  connect(switch1.u1, controller.y) annotation (Line(points={{-26,8},{18,8},{18,0},{36.98,0}},
                                                                                   color={0,0,127}));
  connect(switch1.u3, realExpression9.y) annotation (Line(points={{-26,-8},{-12,-8},{-12,-20},{-11,-20}},     color={0,0,127}));
  connect(switch1.y, abs1.u) annotation (Line(points={{-49,0},{-60,0}},   color={0,0,127}));
  connect(abs1.y, firstOrder1.u) annotation (Line(points={{-83,0},{-90,0}},   color={0,0,127}));
  connect(u_switch, controller.reset) annotation (Line(points={{-100,-136},{-100,-32},{56.02,-32},{56.02,-16.96}},     color={255,0,255}));
  connect(u_switch, switch1.u2) annotation (Line(points={{-100,-136},{-100,-32},{14,-32},{14,0},{-26,0}},    color={255,0,255}));
  connect(firstOrder1.y, y) annotation (Line(points={{-113,0},{-150,0}}, color={0,0,127}));
  connect(controller.u_set, u_s) annotation (Line(points={{73.02,8.32},{108,8.32},{108,0},{144,0}}, color={0,0,127}));
  connect(u_m, controller.u_m) annotation (Line(points={{0,-144},{0,-78},{90,-78},{90,-8.96},{73.7,-8.96}}, color={0,0,127}));
  annotation (Diagram(coordinateSystem(extent={{-140,-140},{140,140}})), Icon(coordinateSystem(extent={{-140,-140},{140,140}}), graphics={Rectangle(
          extent={{-140,140},{140,-140}},
          lineColor={28,108,200},
          lineThickness=0.5,
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid), Text(
          extent={{-102,86},{44,-18}},
          textColor={0,0,0},
          textString="Pump"),             Text(
          extent={{-102,38},{110,-130}},
          textColor={0,0,0},
          textString="Controller")}));
end PumpController;
