within ThermofluidStream.HeatExchangersPhysical.ConductionElements;
model CrossIcon "Icon for ConductionElement Cross"

  annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
        Line(
          points={{-80,0},{84,0}},
          color={28,108,200},
          thickness=0.5),
        Ellipse(
          extent={{56,12},{64,-12}},
          lineColor={28,108,200},
          lineThickness=0.5,
          fillColor={170,213,255},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-60,12},{60,-12}},
          lineColor={238,46,47},
          lineThickness=0.5,
          fillColor={170,213,255},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-60,40},{60,-40}},
          lineColor={0,0,0},
          lineThickness=0.5),
        Line(
          points={{-40,30},{-40,-30}},
          color={238,46,47},
          thickness=1),
        Line(
          points={{40,30},{40,-30}},
          color={238,46,47},
          thickness=1),
        Line(
          points={{20,30},{20,-30}},
          color={238,46,47},
          thickness=1),
        Line(
          points={{0,30},{0,-30}},
          color={238,46,47},
          thickness=1),
        Line(
          points={{-20,30},{-20,-30}},
          color={238,46,47},
          thickness=1),
        Ellipse(
          extent={{-64,12},{-56,-12}},
          lineColor={28,108,200},
          lineThickness=0.5,
          fillColor={170,213,255},
          fillPattern=FillPattern.Solid),
        Line(
          points={{0,88},{0,40}},
          color={238,46,47},
          thickness=0.5),
        Text(
          extent={{-60,-46},{60,-74}},
          textColor={0,0,0},
          textString="Cross Flow")}),
    Diagram(coordinateSystem(preserveAspectRatio=false)),
    Documentation(info="<html>
<p>Implementation of the Conduction Element for the DiscritizedHex.</p>
<p>Concerning the heat transfer it is assumed, that the main term influencing the coefficient of heat transfer is the mass flow rate. Therefore a nominal value for the heat transfer coefficient at a nominal mass flow rate can be set. Furthermore a minimum value U_min for the coefficent of heat transfer is set to ensure heat transfer at zero mass flow.</p>
<p>For further documentation see the documentation of the motherclass.</p>
</html>"));
end CrossIcon;
