within ThermofluidStream.HeatExchangersPhysical.HXcellPairs;
model HXcellPair_CondensingAirCross_2PhasePipe "Heat Exchange in Cell: Condensing Air in Cross Element and 2 Phase Fluid in Pipe Element"

  import Modelica.Units.SI;

  replaceable package Medium_CrossChannel =
      ThermofluidStream.Media.myMedia.Air.MoistAir "Medium model in cross channel (A)"
    annotation (choicesAllMatching=true, Dialog(group = "Medium definitions"));

  replaceable package Medium_Water =
      ThermofluidStream.Media.myMedia.Water.StandardWater "Medium model for condensed vapor in cross channel (A)"
    annotation (choicesAllMatching=true, Dialog(group = "Medium definitions"));

  replaceable package Medium_PipeChannel = ThermofluidStream.Media.myMedia.Interfaces.PartialTwoPhaseMedium "Medium model in pipe channel (B)"
    annotation (choicesAllMatching=true, Dialog(group = "Medium definitions"));                             // e.g. Refrigerant

//  constant SI.Density TinDensity = 8920;         // Density of copper
//  constant SI.SpecificHeatCapacity TinCp = 385;  // cp of copper

  parameter HXutilities.MaterialProperties.Material HX_material=HXutilities.MaterialProperties.Material.Copper;

  parameter ThermofluidStream.Processes.Internal.InitializationMethodsCondElement init_pipe=ThermofluidStream.Processes.Internal.InitializationMethodsCondElement.inlet "Initialization method for h" annotation (Dialog(tab="Initialization", group="Pipe"));
  parameter SI.Temperature T_0_pipe=313.15      "Initial Temperature"
    annotation(Dialog(tab="Initialization", group="Pipe", enable=(init == ThermofluidStream.Processes.Internal.InitializationMethodsCondElement.T)));
  parameter SI.SpecificEnthalpy h_0_pipe= 30000 "Initial specific enthalpy" annotation (Dialog(
      tab="Initialization",
      group="Pipe",
      enable=(init_pipe == ThermofluidStream.Processes.Internal.InitializationMethodsCondElement.h)));

  parameter SI.Temperature T_0_cross             "Initial Temperature"
    annotation(Dialog(tab="Initialization", group="Cross"));

  parameter SI.Length PipeLength = 0.28;
  parameter SI.Length PipeInnerDiameter = 0.007;
  parameter SI.Length PipeWallThickness = 0.0025;
  parameter SI.Length CrossHeight = 0.035;
  parameter Integer   PlateFinNumber = 25;
  parameter SI.Length PlateFinLength = 0.035;
  parameter SI.Length PlateFinHeight = 0.035;
  parameter SI.Length PlateFinThickness = 0.0005;

  parameter SI.Temperature T_pipe_initial=300      "Initial Pipe Channel Temperature";
  parameter SI.Temperature T_cross_initial=300     "Initial Cross Channel Temperature";

  final parameter SI.Density TinDensity =         HXutilities.MaterialProperties.MaterialDensity(
                                                                            HX_material);
  final parameter SI.SpecificHeatCapacity TinCp = HXutilities.MaterialProperties.MaterialSpecHeatCapacity(
                                                                                     HX_material);

  final parameter SI.Length PipeOuterDiameter = PipeInnerDiameter + 2*PipeWallThickness;
  final parameter SI.Area   PipeFlowArea = Modelica.Constants.pi * (PipeInnerDiameter/2)^2;
  final parameter SI.Area   PipeOuterCSArea = Modelica.Constants.pi * (PipeOuterDiameter/2)^2;
  final parameter SI.Area   PipeOuterSurface = Modelica.Constants.pi * PipeOuterDiameter * PipeLength;
  final parameter SI.Volume PipeInnerVolume = PipeFlowArea * PipeLength;

  final parameter SI.Area   PlateFinArea = (PlateFinLength * PlateFinHeight - PipeOuterCSArea) * PlateFinNumber;
  final parameter SI.Area   PlateFinHXArea = 2 * PlateFinArea;

  final parameter SI.Length CrossFlowLength = PlateFinLength;
  final parameter SI.Length CrossHydDiameter = 2*(CrossHeight-PipeOuterDiameter)*PipeLength / (CrossHeight-PipeOuterDiameter+PipeLength);
  final parameter SI.Area   CrossFlowArea = Modelica.Constants.pi * (CrossHydDiameter/2)^2;
  final parameter SI.Area   CrossHXArea = PipeOuterSurface + PlateFinHXArea;
  final parameter SI.Volume CrossFlowVolume = (CrossFlowLength * CrossHeight - PipeOuterCSArea) * PipeLength;

  final parameter SI.Volume PipeWallVolume = (PipeOuterCSArea - PipeFlowArea) * PipeLength;
  final parameter SI.Mass   PipeWallMass = TinDensity * PipeWallVolume;
  final parameter SI.Volume PlateFinVolume = PlateFinArea * PlateFinThickness;
  final parameter SI.Mass   PlateFinMass = TinDensity * PlateFinVolume;
  final parameter SI.Mass   TinMass = PipeWallMass + PlateFinMass;

  final parameter SI.HeatCapacity TinHeatCapacity = TinCp * TinMass;

  SI.MassFlowRate m_flow_pipe_cond = conductionElement_pipe.m_flow_cond;

  SI.HeatFlowRate Q_flow_cross = conductionElement_cross.heatPort.Q_flow;
  SI.HeatFlowRate Q_flow_pipe =  conductionElement_pipe.heatPort.Q_flow;

  ThermofluidStream.Interfaces.Inlet inlet_Pipe(redeclare package Medium = Medium_PipeChannel)
    annotation (Placement(transformation(extent={{-10,-10},{10,10}},
        rotation=90,
        origin={0,-100}), iconTransformation(
        extent={{-10,-10},{10,10}},
        rotation=90,
        origin={0,-100})));
  ThermofluidStream.Interfaces.Outlet outlet_Pipe(redeclare package Medium = Medium_PipeChannel)
    annotation (Placement(transformation(extent={{-10,-10},{10,10}},
        rotation=90,
        origin={0,100}), iconTransformation(
        extent={{-10,-10},{10,10}},
        rotation=90,
        origin={0,100})));
  ThermofluidStream.Interfaces.Inlet inlet_Cross(redeclare package Medium = Medium_CrossChannel) annotation (
      Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={-100,0}), iconTransformation(extent={{-10,-10},{10,10}}, origin={-100,0})));
  ThermofluidStream.Interfaces.Outlet outlet_Cross(redeclare package Medium = Medium_CrossChannel) annotation (
      Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={100,0}), iconTransformation(extent={{-10,-10},{10,10}}, origin={100,0})));
  ConductionElements.ConductionElement_Pipe_2Phase conductionElement_pipe(
    init=ThermofluidStream.Processes.Internal.InitializationMethodsCondElement.T,
    T_0=T_pipe_initial,
    redeclare package Medium = Medium_PipeChannel,
    V=PipeInnerVolume,
    d_hyd=PipeInnerDiameter,
    L_flow=PipeLength) annotation (Placement(transformation(
        extent={{-20,-20},{20,20}},
        rotation=90,
        origin={40,40})));
  ConductionElements.ConductionElement_Cross_CondensingAir conductionElement_cross(
    init=ThermofluidStream.Processes.Internal.InitializationMethodsCondElement.T,
    redeclare package Medium_Water = Medium_Water,
    d_out=PipeOuterDiameter,
    L_pipe=PipeLength,
    H_ch=CrossHeight,
    L_ch=PlateFinLength,
    A_HX=CrossHXArea,
    N_fin=PlateFinNumber,
    redeclare package Medium = Medium_CrossChannel,
    V=CrossFlowVolume,
    T_0=T_0_cross) annotation (Placement(transformation(
        extent={{-20,-20},{20,20}},
        rotation=0,
        origin={-40,-28})));
  Modelica.Thermal.HeatTransfer.Components.HeatCapacitor Tin(C=TinHeatCapacity)
    annotation (Placement(transformation(extent={{-20,-20},{20,20}},
        rotation=0,
        origin={-40,60})));
  ThermofluidStream.Interfaces.Inlet inlet_cross_water(redeclare package Medium = Medium_Water) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={-100,-40}), iconTransformation(extent={{-10,-10},{10,10}}, origin={-100,-40})));
  ThermofluidStream.Interfaces.Outlet outlet_cross_water(redeclare package Medium = Medium_Water) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={100,-40}), iconTransformation(extent={{-10,-10},{10,10}}, origin={100,-40})));
  Modelica.Thermal.HeatTransfer.Sensors.HeatFlowSensor heatFlowSensor annotation (Placement(transformation(extent={{0,30},{-20,50}})));
  Modelica.Thermal.HeatTransfer.Sensors.HeatFlowSensor heatFlowSensor1 annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=90,
        origin={-40,10})));
  Modelica.Thermal.HeatTransfer.Celsius.TemperatureSensor temperatureSensor annotation (Placement(transformation(extent={{-60,30},{-80,50}})));
  Modelica.Blocks.Interaction.Show.RealValue realValue(significantDigits=3) annotation (Placement(transformation(extent={{-20,0},{0,20}})));
  Modelica.Blocks.Interaction.Show.RealValue realValue1(significantDigits=3) annotation (Placement(transformation(extent={{2,14},{22,34}})));
  Modelica.Blocks.Interaction.Show.RealValue realValue2 annotation (Placement(transformation(extent={{-94,30},{-114,50}})));
equation

  connect(conductionElement_pipe.outlet, outlet_Pipe)
    annotation (Line(
      points={{40,60},{40,80},{0,80},{0,100}},
      color={28,108,200},
      thickness=0.5));
  connect(inlet_Cross, conductionElement_cross.inlet) annotation (Line(
      points={{-100,0},{-70,0},{-70,-28},{-60,-28}},
      color={28,108,200},
      thickness=0.5));
  connect(conductionElement_cross.outlet, outlet_Cross) annotation (Line(
      points={{-20,-28},{74,-28},{74,0},{100,0}},
      color={28,108,200},
      thickness=0.5));
  connect(inlet_Pipe, conductionElement_pipe.inlet) annotation (Line(
      points={{0,-100},{0,-80},{40,-80},{40,20}},
      color={28,108,200},
      thickness=0.5));
  connect(inlet_cross_water, conductionElement_cross.inlet_water) annotation (Line(
      points={{-100,-40},{-76,-40},{-76,-40},{-60,-40}},
      color={28,108,200},
      thickness=0.5));
  connect(conductionElement_cross.outlet_water, outlet_cross_water) annotation (Line(
      points={{-20,-40},{46,-40},{46,-40},{100,-40}},
      color={28,108,200},
      thickness=0.5));
  connect(conductionElement_cross.heatPort, heatFlowSensor1.port_a) annotation (Line(points={{-40,-8.4},{-40,0}}, color={191,0,0}));
  connect(heatFlowSensor1.port_b, Tin.port) annotation (Line(points={{-40,20},{-40,40}}, color={191,0,0}));
  connect(temperatureSensor.port, Tin.port) annotation (Line(points={{-60,40},{-40,40}}, color={191,0,0}));
  connect(realValue.numberPort, heatFlowSensor1.Q_flow) annotation (Line(points={{-21.5,10},{-29,10}}, color={0,0,127}));
  connect(heatFlowSensor.Q_flow, realValue1.numberPort) annotation (Line(points={{-10,29},{-10,24},{0.5,24}}, color={0,0,127}));
  connect(realValue2.numberPort, temperatureSensor.T) annotation (Line(points={{-92.5,40},{-80,40}}, color={0,0,127}));
  connect(Tin.port, heatFlowSensor.port_b) annotation (Line(points={{-40,40},{-20,40}}, color={191,0,0}));
  connect(heatFlowSensor.port_a, conductionElement_pipe.heatPort) annotation (Line(points={{0,40},{20.4,40}}, color={191,0,0}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=
            false), graphics={Text(
          extent={{56,94},{96,72}},
          textColor={0,0,0},
          textString="Pipe Flow"), Text(
          extent={{58,-54},{98,-76}},
          textColor={0,0,0},
          textString="Cross Flow")}),
    Documentation(info="<html>
    A general description of cross flow heat exchagers is given in the package.
    
    <p>
    <strong>Cross Flow Media</strong>
    <p>
    Moist Air: condensing or non-condensing<br>
    Condensing Air leaves the cross flow channel in 2 flows: moist air and liquid water
    
    <p>
    <strong>Pipe Flow Media</strong>
    <p>
    1 substance in 1 or 2 phases like refrigerents, water, glycol, thermal oil<br>
    
    </html>"));
end HXcellPair_CondensingAirCross_2PhasePipe;
