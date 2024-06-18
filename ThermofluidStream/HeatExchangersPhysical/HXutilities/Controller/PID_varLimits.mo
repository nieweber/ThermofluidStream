within ThermofluidStream.HeatExchangersPhysical.HXutilities.Controller;
block PID_varLimits
  "P, I, PI, PD, and PID controller with optional variable limits, anti-windup compensation and setpoint weighting"

/* -------------------------------------------------------------- */
/* Copyright                                                      */
/* (c) 2022-2023 DLR Institute of System Dynamics and Control     */
/* All rights reserved: Copying, Publishing, Changing, Exploiting */
/*                                                                */
/* Author(s): Peter Eschenbacher                                  */
/*                                                                */
/* Disclaimer: No Warranty of Any Kind of Errors                  */
/* -------------------------------------------------------------- */

  parameter Controller.ControllerType controllerType=Controller.ControllerType.PI "Type of controller";

  parameter Boolean useMaxPort = false "enable port for maximum output value" annotation(Dialog(group="Limiter"));

  parameter Boolean useMinPort = false "enable port for minimum output value" annotation(Dialog(group="Limiter"));

  parameter Real yMax = 1e10 "maximum output value" annotation(Dialog(group="Limiter", enable=not useMaxPort));

  parameter Real yMin = -1e10 "minimum output value" annotation(Dialog(group="Limiter", enable=not useMinPort));

  parameter Real k(unit="1")=1          "Gain of controller"
    annotation (Dialog(enable=controllerType == Controller.ControllerType.P or controllerType == Controller.ControllerType.PI or controllerType == Controller.ControllerType.PD or controllerType == Controller.ControllerType.PID));
  parameter Modelica.Units.SI.Time Ti(min=Modelica.Constants.small)=0.5 "Time constant of Integrator block"
    annotation (Dialog(enable=controllerType == Controller.ControllerType.I or controllerType == Controller.ControllerType.PI or controllerType == Controller.ControllerType.PID));
  parameter Modelica.Units.SI.Time Td(min=0)=0.1 "Time constant of Derivative block" annotation (Dialog(enable=controllerType == Controller.ControllerType.PD or controllerType == Controller.ControllerType.PID));
  parameter Real wp(min=0) = 1 "Set-point weight for Proportional block (0..1)"
    annotation (Dialog(enable=controllerType == Controller.ControllerType.P or controllerType == Controller.ControllerType.PI or controllerType == Controller.ControllerType.PD or controllerType == Controller.ControllerType.PID));
  parameter Real wd(min=0) = 1 "Set-point weight for Derivative block (0..1)"
     annotation(Dialog(enable=controllerType == Controller.ControllerType.PD or controllerType == Controller.ControllerType.PID));
  parameter Real Ni(min=100*Modelica.Constants.eps) = 0.9
    "Ni*Ti is time constant of anti-windup compensation"
     annotation(Dialog(enable=controllerType == Controller.ControllerType.I or controllerType == Controller.ControllerType.PI or controllerType == Controller.ControllerType.PID));
  parameter Real Nd(min=100*Modelica.Constants.eps) = 10
    "The higher Nd, the more ideal the derivative block"
     annotation(Dialog(enable=controllerType == Controller.ControllerType.PD or controllerType == Controller.ControllerType.PID));
  parameter Real xi_start=0
    "Initial or guess value for integrator output (= integrator state)"
    annotation (Dialog(group="Initialization",
                enable=controllerType == Controller.ControllerType.I or controllerType == Controller.ControllerType.PI or controllerType == Controller.ControllerType.PID));
  parameter Real xd_start=0
    "Initial or guess value for state of derivative block"
    annotation (Dialog(group="Initialization",
                         enable=controllerType == Controller.ControllerType.PD or controllerType == Controller.ControllerType.PID));
  parameter Boolean strict=false "= true, if strict limits with noEvent(..)"
    annotation (Evaluate=true, choices(checkBox=true), Dialog(tab="Advanced"));
  constant Modelica.Units.SI.Time unitTime=1 annotation (HideResult=true);
  Modelica.Blocks.Math.Add addP(k1=wp, k2=-1) if with_P
    annotation (Placement(transformation(extent={{-80,40},{-60,60}})));
  Modelica.Blocks.Math.Add addD(k1=wd, k2=-1) if with_D
    annotation (Placement(transformation(extent={{-80,-10},{-60,10}})));
  Modelica.Blocks.Math.Gain P(k=kP)
                                   if with_P
    annotation (Placement(transformation(extent={{-50,40},{-30,60}})));
  Modelica.Blocks.Continuous.Integrator I(
    k=kI*unitTime/Ti,
    y_start=xi_start,
    initType=Modelica.Blocks.Types.Init.InitialState)
                                                  if with_I
    annotation (Placement(transformation(extent={{-50,-60},{-30,-40}})));
  Modelica.Blocks.Continuous.Derivative D(
    k=kD*Td/unitTime,
    T=max([Td/Nd,1.e-14]),
    x_start=xd_start,
    initType=Modelica.Blocks.Types.Init.InitialState)
                                            if with_D
    annotation (Placement(transformation(extent={{-50,-10},{-30,10}})));
  Modelica.Blocks.Math.Gain gainPID(k=1)
    annotation (Placement(transformation(extent={{20,-10},{40,10}})));
  Modelica.Blocks.Math.Add3 addPID
    annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
  Modelica.Blocks.Math.Add3 addI(k2=-1) if with_I
    annotation (Placement(transformation(extent={{-80,-60},{-60,-40}})));
  Modelica.Blocks.Math.Add addSat(k1=+1, k2=-1) if with_I annotation (Placement(
        transformation(
        origin={70,-50},
        extent={{-10,-10},{10,10}},
        rotation=270)));
  Modelica.Blocks.Math.Gain gainTrack(k=1/(k*Ni)) if with_I
    annotation (Placement(transformation(extent={{0,-80},{-20,-60}})));
protected
  parameter Boolean with_P = controllerType==Controller.ControllerType.P
                                                              or
                             controllerType==Controller.ControllerType.PI
                                                               or
                             controllerType==Controller.ControllerType.PD
                                                               or
                             controllerType==Controller.ControllerType.PID
                                                                annotation(Evaluate=true, HideResult=true);
  parameter Boolean with_I = controllerType==Controller.ControllerType.I
                                                              or
                             controllerType==Controller.ControllerType.PI
                                                               or
                             controllerType==Controller.ControllerType.PID
                                                                annotation(Evaluate=true, HideResult=true);
  parameter Boolean with_D = controllerType==Controller.ControllerType.PD
                                                               or
                             controllerType==Controller.ControllerType.PID
                                                                annotation(Evaluate=true, HideResult=true);

  parameter Real kP = k;
  parameter Real kI = if with_P then k else 1;
  parameter Real kD = k;

public
  Modelica.Blocks.Sources.Constant Dzero(k=0) if not with_D
    annotation (Placement(transformation(extent={{-40,20},{-30,30}})));
  Modelica.Blocks.Sources.Constant Izero(k=0) if not with_I
    annotation (Placement(transformation(extent={{0,-55},{-10,-45}})));
  Modelica.Blocks.Math.Min min annotation (Placement(transformation(extent={{60,-10},{80,10}})));
  Modelica.Blocks.Interfaces.RealInput y_max if useMaxPort
                                             "Connector of maximum output signal" annotation (Placement(
        transformation(
        extent={{-10,10},{10,-10}},
        rotation=180,
        origin={110,80}), iconTransformation(
        extent={{-10,10},{10,-10}},
        rotation=180,
        origin={110,80})));
  Modelica.Blocks.Math.Max max1 annotation (Placement(transformation(extent={{72,42},{92,62}})));
public
  Modelica.Blocks.Sources.Constant Pzero(k=0) if not with_P
    annotation (Placement(transformation(extent={{-40,70},{-30,80}})));
  Modelica.Blocks.Interfaces.RealOutput
             y "Connector of actuator output signal" annotation (Placement(
        transformation(extent={{100,-10},{120,10}})));
  Modelica.Blocks.Interfaces.RealInput
            u_m "Connector of measurement input signal" annotation (Placement(
        transformation(
        origin={0,-110},
        extent={{10,-10},{-10,10}},
        rotation=270), iconTransformation(
        extent={{10,-10},{-10,10}},
        rotation=270,
        origin={0,-110})));
  Modelica.Blocks.Interfaces.RealInput
            u_s "Connector of setpoint input signal" annotation (Placement(
        transformation(extent={{-120,-10},{-100,10}}), iconTransformation(extent={{-120,-10},{-100,10}})));
  Modelica.Blocks.Sources.Constant fixedMaxVal(k=yMax) if not useMaxPort
    annotation (Placement(transformation(extent={{0,80},{20,100}})));
  Modelica.Blocks.Sources.Constant fixedMinVal(k=yMin) if not useMinPort
    annotation (Placement(transformation(extent={{0,48},{20,68}})));
  Modelica.Blocks.Interfaces.RealInput y_min if useMinPort annotation (Placement(transformation(
        extent={{-10,10},{10,-10}},
        rotation=180,
        origin={110,-80})));
initial equation

equation

  connect(u_s, addP.u1) annotation (Line(points={{-110,0},{-94,0},{-94,56},{-82,56}},
                    color={0,0,127}));
  connect(u_s, addD.u1) annotation (Line(points={{-110,0},{-94,0},{-94,6},{-82,6}},
                   color={0,0,127}));
  connect(u_s, addI.u1) annotation (Line(points={{-110,0},{-90,0},{-90,-42},{-82,-42}},
                     color={0,0,127}));
  connect(addP.y, P.u) annotation (Line(points={{-59,50},{-52,50}}, color={0,
          0,127}));
  connect(addD.y, D.u)
    annotation (Line(points={{-59,0},{-52,0}}, color={0,0,127}));
  connect(addI.y, I.u) annotation (Line(points={{-59,-50},{-52,-50}}, color={
          0,0,127}));
  connect(P.y, addPID.u1) annotation (Line(points={{-29,50},{-20,50},{-20,8},{-12,
          8}},     color={0,0,127}));
  connect(D.y, addPID.u2)
    annotation (Line(points={{-29,0},{-12,0}},color={0,0,127}));
  connect(I.y, addPID.u3) annotation (Line(points={{-29,-50},{-20,-50},{-20,-8},
          {-12,-8}},    color={0,0,127}));
  connect(addSat.y, gainTrack.u) annotation (Line(points={{70,-61},{70,-70},{2,-70}},
                    color={0,0,127}));
  connect(gainTrack.y, addI.u3) annotation (Line(points={{-21,-70},{-88,-70},{-88,
          -58},{-82,-58}},     color={0,0,127}));
  connect(u_m, addP.u2) annotation (Line(points={{0,-110},{0,-92},{-92,-92},{-92,44},{-82,44}}, color={0,0,127}));
  connect(u_m, addD.u2) annotation (Line(points={{0,-110},{0,-92},{-92,-92},{-92,-6},{-82,-6}}, color={0,0,127}));
  connect(u_m, addI.u2) annotation (Line(points={{0,-110},{0,-92},{-92,-92},{-92,-50},{-82,-50}}, color={0,0,127}));
  connect(Dzero.y, addPID.u2) annotation (Line(points={{-29.5,25},{-24,25},{-24,
          0},{-12,0}},    color={0,0,127}));
  connect(Izero.y, addPID.u3) annotation (Line(points={{-10.5,-50},{-20,-50},{-20,
          -8},{-12,-8}},    color={0,0,127}));
  connect(addPID.y, gainPID.u)
    annotation (Line(points={{11,0},{18,0}}, color={0,0,127}));
  connect(addSat.u2, gainPID.y)
    annotation (Line(points={{64,-38},{64,-30},{46,-30},{46,0},{41,0}}, color={0,0,127}));
  connect(y_max, min.u1) annotation (Line(points={{110,80},{54,80},{54,6},{58,6}},       color={0,0,127}));
  connect(min.y, addSat.u1) annotation (Line(points={{81,0},{90,0},{90,-30},{76,-30},{76,-38}}, color={0,0,127}));
  connect(min.u2, max1.y)
    annotation (Line(points={{58,-6},{50,-6},{50,36},{100,36},{100,52},{93,52}}, color={0,0,127}));
  connect(max1.u2, gainPID.y) annotation (Line(points={{70,46},{46,46},{46,0},{41,0}}, color={0,0,127}));
  connect(Pzero.y, addPID.u1) annotation (Line(points={{-29.5,75},{-20,75},{-20,8},{-12,8}}, color={0,0,127}));
  connect(min.y, y) annotation (Line(points={{81,0},{110,0}}, color={0,0,127}));
  connect(fixedMaxVal.y, min.u1) annotation (Line(points={{21,90},{54,90},{54,6},{58,6}}, color={0,0,127}));
  connect(fixedMinVal.y, max1.u1) annotation (Line(points={{21,58},{70,58}}, color={0,0,127}));
  connect(y_min, max1.u1) annotation (Line(points={{110,-80},{44,-80},{44,58},{70,58}}, color={0,0,127}));
  annotation (defaultComponentName="PID",
    Icon(coordinateSystem(
        preserveAspectRatio=true,
        extent={{-100,-100},{100,100}}), graphics={
        Line(points={{-80,78},{-80,-90}}, color={192,192,192}),
        Polygon(
          points={{-80,90},{-88,68},{-72,68},{-80,90}},
          lineColor={192,192,192},
          fillColor={192,192,192},
          fillPattern=FillPattern.Solid),
        Line(points={{-90,-80},{82,-80}}, color={192,192,192}),
        Polygon(
          points={{90,-80},{68,-72},{68,-88},{90,-80}},
          lineColor={192,192,192},
          fillColor={192,192,192},
          fillPattern=FillPattern.Solid),
        Line(points={{-80,-80},{-80,-20},{30,60},{80,60}}, color={0,0,127}),
        Text(
          extent={{-20,-20},{80,-60}},
          textColor={192,192,192},
          textString="%controllerType"),
        Line(
          visible=strict,
          points={{30,60},{81,60}},
          color={255,0,0}),
        Rectangle(extent={{-100,100},{100,-100}}, lineColor={28,108,200}),
                                        Text(
        extent={{-150,140},{150,100}},
        textString="%name",
        textColor={0,0,255})}),
    Documentation(info="<html>
<p>
Via parameter <strong>controllerType</strong> either <strong>P</strong>, <strong>PI</strong>, <strong>PD</strong>,
or <strong>PID</strong> can be selected. If, e.g., PI is selected, all components belonging to the
D-part are removed from the block (via conditional declarations).
The example model
<a href=\"modelica://Modelica.Blocks.Examples.PID_Controller\">Modelica.Blocks.Examples.PID_Controller</a>
demonstrates the usage of this controller.
Several practical aspects of PID controller design are incorporated
according to chapter 3 of the book:
</p>

<dl>
<dt>&Aring;str&ouml;m K.J., and H&auml;gglund T.:</dt>
<dd> <strong>PID Controllers: Theory, Design, and Tuning</strong>.
     Instrument Society of America, 2nd edition, 1995.
</dd>
</dl>

<p>
Besides the additive <strong>proportional, integral</strong> and <strong>derivative</strong>
part of this controller, the following features are present:
</p>
<ul>
<li> The output of this controller is limited. If the controller is
     in its limits, anti-windup compensation is activated to drive
     the integrator state to zero.</li>
<li> The high-frequency gain of the derivative part is limited
     to avoid excessive amplification of measurement noise.</li>
<li> Setpoint weighting is present, which allows to weight
     the setpoint in the proportional and the derivative part
     independently from the measurement. The controller will respond
     to load disturbances and measurement noise independently of this setting
     (parameters wp, wd). However, setpoint changes will depend on this
     setting. For example, it is useful to set the setpoint weight wd
     for the derivative part to zero, if steps may occur in the
     setpoint signal.</li>
</ul>

<p>
The parameters of the controller can be manually adjusted by performing
simulations of the closed loop system (= controller + plant connected
together) and using the following strategy:
</p>

<ol>
<li> Set very large limits</li>
<li> Select a <strong>P</strong>-controller and manually enlarge parameter <strong>k</strong>
     (the total gain of the controller) until the closed-loop response
     cannot be improved any more.</li>
<li> Select a <strong>PI</strong>-controller and manually adjust parameters
     <strong>k</strong> and <strong>Ti</strong> (the time constant of the integrator).
     The first value of Ti can be selected, such that it is in the
     order of the time constant of the oscillations occurring with
     the P-controller. If, e.g., vibrations in the order of T=10 ms
     occur in the previous step, start with Ti=0.01 s.</li>
<li> If you want to make the reaction of the control loop faster
     (but probably less robust against disturbances and measurement noise)
     select a <strong>PID</strong>-Controller and manually adjust parameters
     <strong>k</strong>, <strong>Ti</strong>, <strong>Td</strong> (time constant of derivative block).</li>
<li> Set the limit yMin according to your specification.</li>
<li> Perform simulations such that the output of the PID controller
     goes in its limits. Tune <strong>Ni</strong> (Ni*Ti is the time constant of
     the anti-windup compensation) such that the input to the limiter
     block (= limiter.u) goes quickly enough back to its limits.
     If Ni is decreased, this happens faster. If Ni=infinity, the
     anti-windup compensation is switched off and the controller works bad.</li>
</ol>

<p>
<strong>Initialization</strong>
</p>

<p>
The initial valiues of the integrator (I) and differentiator (D) must be initialized.
</p>

</html>"));
end PID_varLimits;
