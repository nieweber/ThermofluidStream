within TFS_Tutorial.Submodels;
model PTC_Heater

  replaceable package Medium = ThermofluidStream.Media.myMedia.Interfaces.PartialMedium "Medium model";

  parameter Real heatingCapacity( unit = "W") = 500 "Heating capacity [W]";

  parameter Modelica.Units.SI.Temperature T_start "Initial temperature";

  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow2
    annotation (Placement(transformation(extent={{-10,-10},{10,10}},
        rotation=270,
        origin={0,-22})));
  Modelica.Blocks.Math.Gain gain(k=heatingCapacity)
                                         annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=270,
        origin={0,20})));
  Modelica.Blocks.Math.BooleanToReal booleanToReal1(realTrue=1, realFalse=0)  annotation (Placement(transformation(extent={{-10,-10},{10,10}},
        rotation=270,
        origin={0,52})));
  Modelica.Thermal.HeatTransfer.Components.HeatCapacitor heatCapacitorBattery1(C=1000, T(start=T_start, fixed=true))
                                                                                                                   annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=90,
        origin={-10,-52})));
  ThermofluidStream.Processes.ConductionElement batteryPTC(
    redeclare package Medium = Medium,
    init=ThermofluidStream.Processes.Internal.InitializationMethodsCondElement.T,
    T_0=283.15,
    A=10) annotation (Placement(transformation(extent={{-10,10},{10,-10}},
        rotation=180,
        origin={0,-80})));
  ThermofluidStream.Interfaces.Inlet inlet(redeclare package Medium = Medium)
                                           annotation (Placement(transformation(extent={{108,-50},{88,-30}}), iconTransformation(extent={{108,-50},{88,-30}})));
  ThermofluidStream.Interfaces.Outlet outlet(redeclare package Medium = Medium)
                                             annotation (Placement(transformation(extent={{-92,-50},{-112,-30}}), iconTransformation(extent={{-92,-50},{-112,-30}})));
  Modelica.Blocks.Interfaces.BooleanInput u annotation (Placement(transformation(
        extent={{-20,-20},{20,20}},
        rotation=270,
        origin={0,106})));
equation
  connect(booleanToReal1.y, gain.u) annotation (Line(points={{-2.22045e-15,41},{2.22045e-15,41},{2.22045e-15,32}},
                                                                                                          color={0,0,127}));
  connect(heatCapacitorBattery1.port, batteryPTC.heatPort) annotation (Line(points={{0,-52},{0,-70.2},{4.44089e-16,-70.2}}, color={191,0,0}));
  connect(prescribedHeatFlow2.port, heatCapacitorBattery1.port) annotation (Line(points={{-1.77636e-15,-32},{0,-32},{0,-52}},
                                                                                                                            color={191,0,0}));
  connect(batteryPTC.inlet, inlet) annotation (Line(
      points={{10,-80},{98,-80},{98,-40}},
      color={28,108,200},
      thickness=0.5));
  connect(batteryPTC.outlet, outlet) annotation (Line(
      points={{-10,-80},{-56,-80},{-56,-40},{-102,-40}},
      color={28,108,200},
      thickness=0.5));
  connect(u, booleanToReal1.u) annotation (Line(points={{0,106},{0,64},{2.22045e-15,64}}, color={255,0,255}));
  connect(gain.y, prescribedHeatFlow2.Q_flow) annotation (Line(points={{0,9},{0,-12}}, color={0,0,127}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
        Rectangle(
          extent={{-100,100},{100,-100}},
          lineColor={0,0,0},
          lineThickness=0.5,
          fillColor=DynamicSelect({255,255,255}, if u then {255,255,255} else {238,46,47}),
          fillPattern=FillPattern.Solid),
       Ellipse(
         extent={{-24,-64},{24,-16}},
         lineColor={28,108,200},
         lineThickness=0.5,
         fillColor={170,213,255},
         fillPattern=FillPattern.Solid,
         pattern=LinePattern.Solid),
       Line(
         points={{0,-2},{0,-50}},
         color={238,46,47}),
       Line(
         points={{-10,-54},{12,-54}},
         color={238,46,47}),
       Line(
         points={{-10,-50},{12,-50}},
         color={238,46,47}),
       Line(
         points={{-10,-42},{12,-42}},
         color={238,46,47}),
       Line(
         points={{-10,-46},{12,-46}},
         color={238,46,47}),
       Line(
         points={{-10,-34},{12,-34}},
         color={238,46,47}),
       Line(
         points={{-10,-38},{12,-38}},
         color={238,46,47}),
        Line(
          points={{90,-40},{24,-40}},
          color={28,108,200},
          thickness=0.5),
        Line(
          points={{-24,-40},{-96,-40}},
          color={28,108,200},
          thickness=0.5),
        Text(
          extent={{-74,80},{70,-4}},
          textColor={0,0,0},
          textString="HEATER")}),                                Diagram(coordinateSystem(preserveAspectRatio=false)));
end PTC_Heater;
