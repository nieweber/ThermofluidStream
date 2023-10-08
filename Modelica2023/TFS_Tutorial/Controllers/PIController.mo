within TFS_Tutorial.Controllers;
model PIController

  parameter Real Ni(min=100*Modelica.Constants.eps) = 0.9
    "Ni*Ti is time constant of anti-windup compensation";

  parameter Real k(unit="1")=1   "Gain of controller";
  parameter Modelica.Units.SI.Time Ti(min=Modelica.Constants.small)=0.5 "Time constant of Integrator block";

  Modelica.Blocks.Continuous.LimIntegrator limIntegrator(
    k=-1/Ti,
    outMax=1000,                                         use_reset=true) annotation (Placement(transformation(extent={{-22,-30},{-2,-10}})));
  Modelica.Blocks.Math.Gain gain(k=k)
                                 annotation (Placement(transformation(extent={{-22,28},{-2,48}})));
  Modelica.Blocks.Math.Add add annotation (Placement(transformation(extent={{38,-10},{58,10}})));
  Modelica.Blocks.Math.Feedback feedback annotation (Placement(transformation(extent={{-90,-10},{-70,10}})));
  Modelica.Blocks.Interfaces.RealInput u_set annotation (Placement(transformation(extent={{-126,32},{-86,72}})));
  Modelica.Blocks.Interfaces.RealInput u_m annotation (Placement(transformation(extent={{-130,-76},{-90,-36}})));
  Modelica.Blocks.Interfaces.RealOutput y annotation (Placement(transformation(extent={{96,-10},{116,10}})));
  Modelica.Blocks.Interfaces.BooleanInput reset annotation (Placement(transformation(
        extent={{-20,-20},{20,20}},
        rotation=90,
        origin={-6,-106})));
  Modelica.Blocks.Nonlinear.Limiter limiter(uMax=1000) annotation (Placement(transformation(extent={{68,-10},{88,10}})));
  Modelica.Blocks.Math.Add addSat(k1=+1, k2=-1) annotation (Placement(
        transformation(
        origin={78,-46},
        extent={{-10,-10},{10,10}},
        rotation=270)));
  Modelica.Blocks.Math.Gain gainTrack(k=1/(k*Ni))
    annotation (Placement(transformation(extent={{30,-78},{10,-58}})));
  Modelica.Blocks.Math.Add add1 annotation (Placement(transformation(extent={{-52,-30},{-32,-10}})));
equation
  connect(gain.y, add.u1) annotation (Line(points={{-1,38},{16,38},{16,6},{36,6}}, color={0,0,127}));
  connect(limIntegrator.y, add.u2) annotation (Line(points={{-1,-20},{8,-20},{8,-6},{36,-6}}, color={0,0,127}));
  connect(feedback.y, gain.u) annotation (Line(points={{-71,0},{-54,0},{-54,38},{-24,38}}, color={0,0,127}));
  connect(u_set, feedback.u1) annotation (Line(points={{-106,52},{-74,52},{-74,26},{-96,26},{-96,0},{-88,0}}, color={0,0,127}));
  connect(u_m, feedback.u2) annotation (Line(points={{-110,-56},{-80,-56},{-80,-8}}, color={0,0,127}));
  connect(reset, limIntegrator.reset) annotation (Line(points={{-6,-106},{-6,-32}}, color={255,0,255}));
  connect(add.y, limiter.u) annotation (Line(points={{59,0},{66,0}}, color={0,0,127}));
  connect(y, limiter.y) annotation (Line(points={{106,0},{89,0}}, color={0,0,127}));
  connect(addSat.y,gainTrack. u) annotation (Line(points={{78,-57},{78,-68},{32,-68}},
                    color={0,0,127}));
  connect(addSat.u2, add.y) annotation (Line(points={{72,-34},{72,-24},{62,-24},{62,0},{59,0}}, color={0,0,127}));
  connect(addSat.u1, y) annotation (Line(points={{84,-34},{84,-24},{92,-24},{92,0},{106,0}}, color={0,0,127}));
  connect(limIntegrator.u, add1.y) annotation (Line(points={{-24,-20},{-31,-20}}, color={0,0,127}));
  connect(gainTrack.y, add1.u2) annotation (Line(points={{9,-68},{-62,-68},{-62,-26},{-54,-26}}, color={0,0,127}));
  connect(feedback.y, add1.u1) annotation (Line(points={{-71,0},{-62,0},{-62,-14},{-54,-14}}, color={0,0,127}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={Rectangle(
          extent={{-100,100},{100,-100}},
          lineColor={28,108,200},
          lineThickness=0.5,
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid), Text(
          extent={{-78,54},{78,-50}},
          textColor={0,0,0},
          textString="PI")}),                                    Diagram(coordinateSystem(preserveAspectRatio=false)));
end PIController;
