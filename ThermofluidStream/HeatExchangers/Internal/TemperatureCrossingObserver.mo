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

  //DT1 = TinA - ToutB;
  //DT2 = ToutA - TinB;

  DT1 = TinB - ToutA;
  DT2 = ToutB - TinA;

  temperatureCrossing = DT1*DT2;

  annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
          Rectangle(
          extent={{-98,100},{102,-100}},
          lineColor={28,108,200},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid),
        Polygon(points={{-72,-54},{-16,-14},{-10,-20},{-58,-68},{-72,-54}},
                                                                        lineColor={0,0,0},
          fillColor={135,135,135},
          fillPattern=FillPattern.Solid),
        Ellipse(extent={{-28,74},{80,-30}}, lineColor={0,0,0},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{60,-8},{-2,54},{-2,50},{-6,50},{46,-2},{36,-4},{60,-8}},
          pattern=LinePattern.None,
          lineThickness=1,
          fillColor={238,46,47},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{60,54},{-2,-8},{-2,-4},{-6,-4},{46,48},{36,50},{60,54}},
          pattern=LinePattern.None,
          lineThickness=1,
          fillColor={28,108,200},
          fillPattern=FillPattern.Solid)}),
      Diagram(coordinateSystem(preserveAspectRatio=false)));
end TemperatureCrossingObserver;
