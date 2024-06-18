within ThermofluidStream.HeatExchangersPhysical.HXcellPairs;
model HXcellPair_1PhaseCross_CondensingAirPipe "Heat Exchange in Cell: Single Phase Fluid in Cross Element and Condensing Air in Pipe Element"

  import Modelica.Units.SI;

  replaceable package Medium_CrossChannel =
      ThermofluidStream.Media.myMedia.Interfaces.PartialMedium "Medium model in cross channel (A)"
    annotation (choicesAllMatching=true, Dialog(group = "Medium definitions"));

  replaceable package Medium_PipeChannel = ThermofluidStream.Media.myMedia.Air.MoistAir;

  replaceable package Medium_Water =
      ThermofluidStream.Media.myMedia.Water.StandardWater "Medium model for condensed vapor in pipe channel"
    annotation (choicesAllMatching=true, Dialog(group = "Medium definitions"));

//  replaceable package Medium_PipeChannel =
//      ThermofluidStream.Media.myMedia.Interfaces.PartialCondensingGases "Medium model in pipe channel (B)"
//    annotation (choicesAllMatching=true, Dialog(group = "Medium definitions"));                             // Moist Air

  parameter HXutilities.MaterialProperties.Material HX_material=HXutilities.MaterialProperties.Material.Copper;

  parameter SI.Length PipeLength = 0.28;
  parameter SI.Length PipeInnerDiameter = 0.007;
  parameter SI.Length PipeWallThickness = 0.0025;
  parameter SI.Length CrossHeight = 0.035;
  parameter Integer   PlateFinNumber = 25;
  parameter SI.Length PlateFinLength = 0.035;
  parameter SI.Length PlateFinHeight = 0.035;
  parameter SI.Length PlateFinThickness = 0.0005;

  parameter SI.Temperature T_pipe_initial =  300   "Initial Pipe Channel Temperature";
  parameter SI.Temperature T_cross_initial = 300   "Initial Cross Channel Temperature";

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

  SI.MassFlowRate m_flow_pipe_cond = conductionElement_pipe.m_flow_cond_cell;

  SI.MassFlowRate m_flow_pipe_in =        conductionElement_pipe.inlet.m_flow;
  SI.MassFlowRate m_flow_pipe_water_in =  conductionElement_pipe.inlet_water.m_flow;
  SI.MassFlowRate m_flow_pipe_out =       conductionElement_pipe.outlet.m_flow;
  SI.MassFlowRate m_flow_pipe_water_out = conductionElement_pipe.outlet_water.m_flow;
  SI.MassFlowRate m_flow_pipe_balance;

  SI.HeatFlowRate Q_flow_cross = conductionElement_cross.heatPort.Q_flow;
  SI.HeatFlowRate Q_flow_pipe =  conductionElement_pipe.heatPort.Q_flow;

  SI.EnthalpyFlowRate H_flow_cross_in =  conductionElement_cross.H_flow_in;
  SI.EnthalpyFlowRate H_flow_cross_out = conductionElement_cross.H_flow_out;

  SI.EnthalpyFlowRate H_flow_pipe_in =    conductionElement_pipe.H_flow_in;
  SI.EnthalpyFlowRate H_flow_pipe_out =   conductionElement_pipe.H_flow_out;
  SI.EnthalpyFlowRate H_flow_pipe_water = conductionElement_pipe.H_flow_water;
  SI.EnthalpyFlowRate H_flow_balance;

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
  ConductionElements.ConductionElement_Pipe_CondensingAir conductionElement_pipe(
    init=ThermofluidStream.Processes.Internal.InitializationMethodsCondElement.T,
    T_0=T_pipe_initial,
    redeclare package Medium = Medium_PipeChannel,
    V=PipeInnerVolume,
    d_hyd=PipeInnerDiameter,
    L_flow=PipeLength) annotation (Placement(transformation(
        extent={{-20,-20},{20,20}},
        rotation=90,
        origin={0,-50})));
  ConductionElements.ConductionElement_Cross_1Phase conductionElement_cross(
    redeclare package Medium = Medium_CrossChannel,
    V=CrossFlowVolume,
    T_0=T_cross_initial,
    d_out=PipeOuterDiameter,
    L_pipe=PipeLength,
    H_ch=CrossHeight,
    N_fin=PlateFinNumber) annotation (Placement(transformation(
        extent={{-20,20},{20,-20}},
        rotation=0,
        origin={-40,0})));
  Modelica.Thermal.HeatTransfer.Components.HeatCapacitor Tin(C=TinHeatCapacity)
    annotation (Placement(transformation(extent={{-20,20},{20,-20}},
        rotation=0,
        origin={-40,-70})));
  ThermofluidStream.Interfaces.Inlet inlet_pipe_water(redeclare package Medium = Medium_Water) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=90,
        origin={40,-100}), iconTransformation(extent={{-10,-10},{10,10}}, origin={40,-100},
        rotation=90)));
  ThermofluidStream.Interfaces.Outlet outlet_pipe_water(redeclare package Medium = Medium_Water) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=90,
        origin={40,98}),  iconTransformation(extent={{-10,-10},{10,10}}, origin={40,98},
        rotation=90)));
  Modelica.Blocks.Interaction.Show.RealValue H_flow_balance_value(
    use_numberPort=false,
    number=H_flow_balance,
    significantDigits=1) annotation (Placement(transformation(extent={{40,10},{60,30}})));
  Modelica.Blocks.Interaction.Show.RealValue m_flow_pipe_balance_value(
    use_numberPort=false,
    number=m_flow_pipe_balance,
    significantDigits=1) annotation (Placement(transformation(extent={{40,30},{60,50}})));
equation

  m_flow_pipe_balance = (m_flow_pipe_in + m_flow_pipe_water_in) + (m_flow_pipe_out + m_flow_pipe_water_out);

  H_flow_balance = (H_flow_cross_in - H_flow_cross_out) + (H_flow_pipe_in - H_flow_pipe_out - H_flow_pipe_water);

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
  connect(inlet_pipe_water, conductionElement_pipe.inlet_water) annotation (Line(
      points={{40,-100},{40,-86},{12,-86},{12,-70}},
      color={28,108,200},
      thickness=0.5));
  connect(conductionElement_pipe.outlet_water, outlet_pipe_water) annotation (Line(
      points={{12,-30},{12,70},{40,70},{40,98}},
      color={28,108,200},
      thickness=0.5));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=
            false), graphics={Text(
          extent={{-66,90},{-26,68}},
          textColor={0,0,0},
          textString="Pipe Flow"), Text(
          extent={{50,-14},{90,-36}},
          textColor={0,0,0},
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
    Moist Air: condensing or non-condensing<br>
    Condensing Air leaves the pipe in 2 flows: moist air and liquid water
    
    </html>"));
end HXcellPair_1PhaseCross_CondensingAirPipe;
