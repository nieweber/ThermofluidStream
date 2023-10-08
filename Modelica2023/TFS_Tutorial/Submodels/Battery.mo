within TFS_Tutorial.Submodels;
model Battery

  parameter Modelica.Units.SI.Temperature ambientTemperature "Ambient temperature";

  Modelica.Thermal.HeatTransfer.Components.HeatCapacitor heatCapacitorBattery(C=1e5, T(start=ambientTemperature, fixed=true))
                                                                                                                 annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=270,
        origin={10,0})));
  Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermalConductor(G=5) annotation (Placement(transformation(
        extent={{10,-10},{-10,10}},
        rotation=0,
        origin={-44,0})));
  Modelica.Thermal.HeatTransfer.Sources.FixedTemperature      fixedTemperature(T=ambientTemperature)
                                                                                    annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={-76,0})));
  Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a port_a annotation (Placement(transformation(extent={{-10,-110},{10,-90}})));
equation
  connect(heatCapacitorBattery.port,thermalConductor. port_a) annotation (Line(points={{1.77636e-15,1.77636e-15},{1.77636e-15,0},{-34,0}},
                                                                                                                      color={191,0,0}));
  connect(thermalConductor.port_b, fixedTemperature.port) annotation (Line(points={{-54,0},{-66,0}},   color={191,0,0}));
  connect(heatCapacitorBattery.port, port_a) annotation (Line(points={{0,0},{0,-100},{0,-100}}, color={191,0,0}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={Rectangle(
          extent={{-100,100},{102,-100}},
          lineColor={28,108,200},
          lineThickness=0.5,
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid), Text(
          extent={{-54,64},{62,-62}},
          textColor={0,0,0},
          textString="BAT")}),                                   Diagram(coordinateSystem(preserveAspectRatio=false)));
end Battery;
