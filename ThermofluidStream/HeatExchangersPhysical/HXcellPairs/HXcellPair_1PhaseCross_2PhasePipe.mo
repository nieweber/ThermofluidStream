within ThermofluidStream.HeatExchangersPhysical.HXcellPairs;
model HXcellPair_1PhaseCross_2PhasePipe  "Heat Exchange in Cell: Single Phase Fluid in Cross Element and 2 Phase Fluid in Pipe Element"

  import Modelica.Units.SI;

  replaceable package Medium_CrossChannel =
      ThermofluidStream.Media.myMedia.Interfaces.PartialMedium "Medium model in cross channel (A)"
    annotation (choicesAllMatching=true, Dialog(group = "Medium definitions"));
  replaceable package Medium_PipeChannel =
      ThermofluidStream.Media.myMedia.Interfaces.PartialTwoPhaseMedium "Medium model in pipe channel (B)"
    annotation (choicesAllMatching=true, Dialog(group = "Medium definitions"));                             // Moist Air
//  replaceable package Medium_PipeChannel =
//      ThermofluidStream.Media.myMedia.Interfaces.PartialTwoPhaseMedium "Medium model in pipe channel (B)"  // Refrigerant
//    annotation (choicesAllMatching=true, Dialog(group = "Medium definitions"));

//  constant SI.Density TinDensity = 8920;         // Density of copper
//  constant SI.SpecificHeatCapacity TinCp = 385;  // cp of copper

  parameter HXutilities.MaterialProperties.Material HX_material=HXutilities.MaterialProperties.Material.Copper;

  parameter ThermofluidStream.Processes.Internal.InitializationMethodsCondElement init_pipe=ThermofluidStream.Processes.Internal.InitializationMethodsCondElement.inlet "Initialization method for h" annotation (Dialog(tab="Initialization", group="Pipe"));
  parameter SI.Temperature T_0_pipe = 40+273.15 "Initial Temperature"
    annotation(Dialog(tab="Initialization", group="Pipe", enable=(init == ThermofluidStream.Processes.Internal.InitializationMethodsCondElement.T)));
  parameter SI.SpecificEnthalpy h_0_pipe= 30000 "Initial specific enthalpy" annotation (Dialog(
      tab="Initialization",
      group="Pipe",
      enable=(init_pipe == ThermofluidStream.Processes.Internal.InitializationMethodsCondElement.h)));

  parameter SI.Temperature T_0_cross = 20+273.15 "Initial Temperature"
    annotation(Dialog(tab="Initialization", group="Cross"));

  parameter SI.Length PipeLength = 0.28;
  parameter SI.Length PipeInnerDiameter = 0.007;
  parameter SI.Length PipeWallThickness = 0.0025;
  parameter SI.Length CrossHeight = 0.035;
  parameter Integer   PlateFinNumber = 25;
  parameter SI.Length PlateFinLength = 0.035;
  parameter SI.Length PlateFinHeight = 0.035;
  parameter SI.Length PlateFinThickness = 0.0005;

//  parameter SI.SpecificEnthalpy h_pipe_initial =  30000   "Initial Pipe Channel Specific Enthalpy";
//  parameter SI.Temperature      T_cross_initial =   300   "Initial Cross Channel Temperature";

  final parameter SI.Density TinDensity =         HXutilities.MaterialProperties.MaterialDensity(
                                                                            HX_material);
  final parameter SI.SpecificHeatCapacity TinCp = HXutilities.MaterialProperties.MaterialSpecHeatCapacity(
                                                                                     HX_material);

  final parameter SI.Length PipeOuterDiameter = PipeInnerDiameter + 2*PipeWallThickness;
  final parameter SI.Area   PipeFlowArea = Modelica.Constants.pi * (PipeInnerDiameter/2)^2;
  final parameter SI.Area   PipeOuterCSArea = Modelica.Constants.pi * (PipeOuterDiameter/2)^2;
//  final parameter SI.Area   PipeInnerSurface = Modelica.Constants.pi * PipeInnerDiameter * PipeLength;
  final parameter SI.Area   PipeOuterSurface = Modelica.Constants.pi * PipeOuterDiameter * PipeLength;
  final parameter SI.Volume PipeInnerVolume = PipeFlowArea * PipeLength;

  final parameter SI.Area   PlateFinArea = (PlateFinLength * PlateFinHeight - PipeOuterCSArea) * PlateFinNumber;
  final parameter SI.Area   PlateFinHXArea = 2 * PlateFinArea;

  final parameter SI.Length CrossFlowLength = PlateFinLength;
//  final parameter SI.Length CrossHydDiameter = 2*(CrossHeight-PipeOuterDiameter)*PipeLength / (CrossHeight-PipeOuterDiameter+PipeLength);
//  final parameter SI.Area   CrossFlowArea = Modelica.Constants.pi * (CrossHydDiameter/2)^2;
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
    annotation (Placement(transformation(extent={{-20,-20},{20,20}},
        rotation=90,
        origin={0,-100}), iconTransformation(
        extent={{-20,-20},{20,20}},
        rotation=90,
        origin={0,-100})));
  ThermofluidStream.Interfaces.Outlet outlet_Pipe(redeclare package Medium = Medium_PipeChannel)
    annotation (Placement(transformation(extent={{-20,-20},{20,20}},
        rotation=90,
        origin={0,100}), iconTransformation(
        extent={{-20,-20},{20,20}},
        rotation=90,
        origin={0,100})));
  ThermofluidStream.Interfaces.Inlet inlet_Cross(redeclare package Medium = Medium_CrossChannel) annotation (
      Placement(transformation(
        extent={{-20,-20},{20,20}},
        rotation=0,
        origin={-100,0}), iconTransformation(extent={{-20,-20},{20,20}}, origin={-100,0})));
  ThermofluidStream.Interfaces.Outlet outlet_Cross(redeclare package Medium = Medium_CrossChannel) annotation (
      Placement(transformation(
        extent={{-20,-20},{20,20}},
        rotation=0,
        origin={100,0}), iconTransformation(extent={{-20,-20},{20,20}}, origin={100,0})));
  ConductionElements.ConductionElement_Pipe_2Phase conductionElement_pipe(
    initM_flow=ThermofluidStream.Utilities.Types.InitializationMethods.state,
    init=init_pipe,
    T_0=T_0_pipe,
    h_0=h_0_pipe,
    redeclare package Medium = Medium_PipeChannel,
    V=PipeInnerVolume,
    d_hyd=PipeInnerDiameter,
    L_flow=PipeLength) annotation (Placement(transformation(
        extent={{-20,-20},{20,20}},
        rotation=90,
        origin={0,-50})));
  ConductionElements.ConductionElement_Cross_1Phase conductionElement_cross(
    redeclare package Medium = Medium_CrossChannel,
    initM_flow=ThermofluidStream.Utilities.Types.InitializationMethods.state,
    V=CrossFlowVolume,
    init=ThermofluidStream.Processes.Internal.InitializationMethodsCondElement.T,
    T_0=T_0_cross,
    H_ch=CrossHeight,
    N_fin=PlateFinNumber,
    d_out=PipeOuterDiameter,
    L_pipe=PipeLength) annotation (Placement(transformation(
        extent={{-20,20},{20,-20}},
        rotation=0,
        origin={-40,0})));
  Modelica.Thermal.HeatTransfer.Components.HeatCapacitor Tin(C=TinHeatCapacity)
    annotation (Placement(transformation(extent={{-20,20},{20,-20}},
        rotation=0,
        origin={-40,-70})));
equation

  connect(conductionElement_pipe.outlet, outlet_Pipe)
    annotation (Line(
      points={{1.33227e-15,-30},{1.33227e-15,35},{0,35},{0,100}},
      color={28,108,200},
      thickness=0.5));
  connect(inlet_Cross, conductionElement_cross.inlet) annotation (Line(
      points={{-100,0},{-60,0}},
      color={28,108,200},
      thickness=0.5));
  connect(conductionElement_cross.outlet, outlet_Cross) annotation (Line(
      points={{-20,0},{100,0}},
      color={28,108,200},
      thickness=0.5));
  connect(conductionElement_pipe.heatPort, Tin.port)
    annotation (Line(points={{-19.6,-50},{-40,-50}},
                                                   color={191,0,0}));
  connect(Tin.port, conductionElement_cross.heatPort)
    annotation (Line(points={{-40,-50},{-40,-19.6}},
                                                   color={191,0,0}));
  connect(inlet_Pipe, conductionElement_pipe.inlet) annotation (Line(
      points={{0,-100},{0,-70}},
      color={28,108,200},
      thickness=0.5));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=
            false), graphics={Text(
          extent={{-66,90},{-26,68}},
          textColor={28,108,200},
          textString="Pipe Flow"), Text(
          extent={{50,-14},{90,-36}},
          textColor={28,108,200},
          textString="Cross Flow")}),
    Documentation(info="<html> 
    A general description of cross flow heat exchagers is given in the package <b>HXcellPairs</b>.
    
    <p>
    <strong>Cross Flow Media</strong>
    <p>
    substance in 1 phase like liquid water, glycol, thermal oil, non-condensing gas<br>
    
    <p>
    <strong>Pipe Flow Media</strong>
    <p>
    1 substance in 1 or 2 phases like refrigerents, water, glycol, thermal oil<br>
    
    </html>"));
end HXcellPair_1PhaseCross_2PhasePipe;
