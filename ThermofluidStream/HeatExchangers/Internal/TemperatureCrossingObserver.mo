within ThermofluidStream.HeatExchangers.Internal;
model TemperatureCrossingObserver

  Modelica.Blocks.Interfaces.RealOutput temperatureCrossing
    annotation (Placement(transformation(extent={{100,-10},{120,10}})));

  Modelica.Blocks.Interfaces.RealOutput DT1(unit = "K")
    annotation (Placement(transformation(extent={{100,40},{120,60}}), iconTransformation(extent={{100,40},{120,60}})));
  Modelica.Blocks.Interfaces.RealOutput DT2(unit = "K")
    annotation (Placement(transformation(extent={{100,-60},{120,-40}}), iconTransformation(extent={{100,-60},{120,-40}})));

  Modelica.Blocks.Interfaces.RealInput TinA(unit = "K") "Temperature at inlet A"
    annotation (Placement(transformation(extent={{-120,-50},{-80,-10}})));
  Modelica.Blocks.Interfaces.RealInput ToutA(unit = "K") "Temperature at outlet A"
    annotation (Placement(transformation(extent={{-120,-90},{-80,-50}})));
  Modelica.Blocks.Interfaces.RealInput TinB(unit = "K") "Temperature at inlet B"
    annotation (Placement(transformation(extent={{-120,50},{-80,90}})));
  Modelica.Blocks.Interfaces.RealInput ToutB(unit = "K") "Temperature at outlet B"
    annotation (Placement(transformation(extent={{-120,10},{-80,50}})));
equation

  DT1 = TinA - ToutB;
  DT2 = ToutA - TinB;

  temperatureCrossing = DT1*DT2;

  annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
          Rectangle(
          extent={{-100,100},{100,-100}},
          lineColor={28,108,200},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid),
        Ellipse(extent={{-10,46},{60,-24}}, lineColor={28,108,200}),
        Polygon(points={{-50,-36},{-6,-6},{0,-14},{-44,-44},{-50,-36}}, lineColor={28,108,200}),
        Line(points={{12,26},{40,-2}}, color={28,108,200}),
        Line(points={{30,0},{40,-2}}, color={28,108,200}),
        Line(
          points={{-3,-3},{7,-1}},
          color={255,0,0},
          origin={31,29},
          rotation=360),
        Line(
          points={{-14,14},{14,-14}},
          color={255,0,0},
          origin={24,14},
          rotation=90),
        Polygon(
          points={{46,30},{8,-8},{10,-2},{4,-4},{32,24},{22,26},{46,30}},
          pattern=LinePattern.None,
          lineThickness=1,
          fillColor={28,108,200},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{46,-8},{8,30},{10,24},{4,26},{32,-2},{22,-4},{46,-8}},
          pattern=LinePattern.None,
          lineThickness=1,
          fillColor={238,46,47},
          fillPattern=FillPattern.Solid)}),
      Diagram(coordinateSystem(preserveAspectRatio=false)));
end TemperatureCrossingObserver;
