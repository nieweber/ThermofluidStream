within ThermofluidStream.Boundaries;
model TerminalSink "Zero mass flow rate sink"

  extends ThermofluidStream.Utilities.DropOfCommonsPlus;

  replaceable package Medium = Media.myMedia.Interfaces.PartialMedium "Medium model"
    annotation (choicesAllMatching=true, Documentation(info="<html>
<p>
Medium package used in the Source. Make sure it is the same as the one
the inlet the source is connected to.
</p>
</html>"));

  Interfaces.Inlet inlet(redeclare package Medium=Medium)
    annotation (Placement(transformation(extent={{-120,-20},{-80,20}})));

equation
  inlet.m_flow = 0;

  annotation (Icon(coordinateSystem(preserveAspectRatio=true), graphics={
       Text(visible=displayInstanceName,
          extent={{-150,60},{150,100}},
          textString="%name",
          textColor=dropOfCommons.instanceNameColor),
        Line(
          points={{-100,0},{-50,0}},
          color={28,108,200},
          thickness=0.5),
        Rectangle(
          extent={{-56,26},{-16,-34}},
          lineColor={28,108,200},
          lineThickness=0.5,
          fillColor={215,215,215},
          fillPattern=FillPattern.Solid,
          pattern=LinePattern.None),
        Rectangle(
          extent={{-60,30},{-20,-30}},
          lineColor={28,108,200},
          lineThickness=0.5,
          fillColor={170,213,255},
          fillPattern=FillPattern.Solid,
          pattern=LinePattern.None),
        Line(
          points={{-20,30},{-60,-30}},
          color={28,108,200},
          thickness=0.5),
        Line(
          points={{-60,30},{-20,-30}},
          color={28,108,200},
          thickness=0.5)}), Diagram(coordinateSystem(preserveAspectRatio=true)),
    Documentation(info="<html>
<p>Sink that terminates the flow. </p>
<p>It imposes a m_flow=0 boundary.</p>
</html>"));
end TerminalSink;
