within ThermofluidStream.HeatExchangersPhysical.ConductionElements;
model PipeIcon

  annotation (
    Icon(
      coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{100,100}}),
      graphics={
        Ellipse(
          extent={{-78,30},{-38,-30}},
          lineColor={28,108,200},
          lineThickness=0.5,
          fillColor={170,213,255},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-58,30},{60,-30}},
          lineColor={28,108,200},
          lineThickness=0.5,
          fillColor={170,213,255},
          fillPattern=FillPattern.Solid,
          pattern=LinePattern.None),
       Line(
         points={{-58,-30},{56,-30}},
         color={238,46,47},
          thickness=1),
       Line(
         points={{-58,30},{56,30}},
         color={238,46,47},
          thickness=1),
       Line(
         points={{0,104},{0,30}},
         color={238,46,47}),
        Ellipse(
          extent={{40,30},{80,-30}},
          lineColor={28,108,200},
          lineThickness=0.5,
          fillColor={170,213,255},
          fillPattern=FillPattern.Solid),
        Text(
          extent={{-60,-34},{60,-68}},
          textColor={0,0,0},
          textString="Pipe Flow")}),
    Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}})));
end PipeIcon;
