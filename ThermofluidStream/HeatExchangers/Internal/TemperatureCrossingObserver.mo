within ThermofluidStream.HeatExchangers.Internal;
model TemperatureCrossingObserver

  Modelica.Blocks.Interfaces.RealOutput temperatureCrossing
    annotation (Placement(transformation(extent={{100,-10},{120,10}})));

  SI.TemperatureDifference DT1;
  SI.TemperatureDifference DT2;

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
        Line(points={{40,-2},{38,8}}, color={28,108,200}),
        Line(points={{30,0},{40,-2}}, color={28,108,200}),
        Line(
          points={{-3,-3},{7,-1}},
          color={28,108,200},
          origin={31,29},
          rotation=360),
        Line(
          points={{-3,-4},{-1,6}},
          color={28,108,200},
          origin={39,22},
          rotation=360),
        Line(
          points={{-14,14},{14,-14}},
          color={28,108,200},
          origin={24,14},
          rotation=90)}),
      Diagram(coordinateSystem(preserveAspectRatio=false)));
end TemperatureCrossingObserver;
