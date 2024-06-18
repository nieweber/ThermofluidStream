within ThermofluidStream.HeatExchangersPhysical;
model HX_1PhaseCross_2PhasePipe "Heat Exchanger: Single Phase Fluid in Cross Flow and 2 Phase Fluid in Pipe Flow"
  extends HXcrossIcon;

  import Modelica.Units.SI;

  replaceable package Medium_CrossChannel =
      ThermofluidStream.Media.myMedia.Interfaces.PartialMedium "Medium model in cross channel (A)"
    annotation (choicesAllMatching=true, Dialog(group = "Medium definitions"));
  replaceable package Medium_PipeChannel =
      ThermofluidStream.Media.myMedia.Interfaces.PartialTwoPhaseMedium "Medium model in pipe channel (B)"
    annotation (choicesAllMatching=true, Dialog(group = "Medium definitions"));

  parameter HXutilities.MaterialProperties.Material HX_material=HXutilities.MaterialProperties.Material.Copper;

  parameter Integer nRows = 8 "Number of rows"
    annotation (Dialog(group="Geometry"));
  parameter Integer nCellPerRow = 1 "Number of cells per row"
    annotation (Dialog(group="Geometry"));
  parameter Integer nPlateFinsPerRow = 25 "Number of plate fins per row"
    annotation (Dialog(group="Geometry"));
  parameter SI.Length length_HX = 0.035 "length of heat exchanger = length of fins"
    annotation (Dialog(group="Geometry"));
  parameter SI.Length width_HX = 0.28 "width of heat exchanger"
    annotation (Dialog(group="Geometry"));
  parameter SI.Length height_HX = 0.28 "height of heat exchanger"
    annotation (Dialog(group="Geometry"));
  parameter SI.Length inner_diameter_pipe = 0.007 "inner diameter of pipe channel"
    annotation (Dialog(group="Geometry"));
  parameter SI.Length thickness_pipe = 0.002 "thickness of pipe wall"
    annotation (Dialog(group="Geometry"));
  parameter SI.Length thickness_plateFins = 0.0005 "thickness of plate fins"
    annotation (Dialog(group="Geometry"));

//  parameter SI.SpecificEnthalpy h_pipe_initial =  30000   "Initial Pipe Channel Specific Enthalpy";
//  parameter SI.Temperature T_pipe_initial =  300   "Initial Pipe Channel Temperature";
//  parameter SI.Temperature T_cross_initial = 300   "Initial Cross Channel Temperature";

  parameter Boolean initializeMassFlow=false "Initialize mass flow at inlets?" annotation(Dialog(tab = "Initialization", group = "Mass flow"));
  parameter Modelica.Units.SI.MassFlowRate m_flow_0_A=0 "Initial mass flow for side A"
    annotation (Dialog(
      tab="Initialization",
      group="Mass flow",
      enable=initializeMassFlow));
  parameter Modelica.Units.SI.MassFlowRate m_flow_0_B=0 "Initial mass flow for side B"
    annotation (Dialog(
      tab="Initialization",
      group="Mass flow",
      enable=initializeMassFlow));

  parameter Modelica.Units.SI.MassFlowRate m_flow_assert(max=0)=-dropOfCommons.m_flow_reg
    "Assertion threshold for negative massflows" annotation (Dialog(tab="Advanced"));
  parameter Boolean enforce_global_energy_conservation = false "If true, exact global energy conservation is enforced by feeding back all energy stored locally back in the system"
    annotation(Dialog(tab="Advanced"));

  parameter Boolean calculate_efficency= false "Enable calculation of efficency";

  parameter ThermofluidStream.Processes.Internal.InitializationMethodsCondElement init_pipe=ThermofluidStream.Processes.Internal.InitializationMethodsCondElement.inlet "Initialization method for h" annotation (Dialog(tab="Initialization", group="Pipe"));
  parameter SI.Temperature T_0_pipe = 40+273.15 "Initial Temperature"
    annotation(Dialog(tab="Initialization", group="Pipe", enable=(init == ThermofluidStream.Processes.Internal.InitializationMethodsCondElement.T)));
  parameter SI.SpecificEnthalpy h_0_pipe= 30000 "Initial specific enthalpy" annotation (Dialog(
      tab="Initialization",
      group="Pipe",
      enable=(init_pipe == ThermofluidStream.Processes.Internal.InitializationMethodsCondElement.h)));

  parameter SI.Temperature T_0_cross = 20+273.15 "Initial Temperature"
    annotation(Dialog(tab="Initialization", group="Cross"));

  final parameter Integer nCells = nRows*nCellPerRow  "Number of discretization elements";

public
  HXcellPairs.HXcellPair_1PhaseCross_2PhasePipe hxCellPair[nCells](
    redeclare package Medium_CrossChannel = Medium_CrossChannel,
    redeclare package Medium_PipeChannel = Medium_PipeChannel,
    each PipeLength=width_HX/nCellPerRow,
    each PipeInnerDiameter=inner_diameter_pipe,
    each CrossHeight=height_HX/nRows,
    each PipeWallThickness=thickness_pipe,
    each PlateFinNumber=div(nPlateFinsPerRow, nCellPerRow),
    each PlateFinLength=length_HX,
    each PlateFinHeight=height_HX/nRows,
    each PlateFinThickness=thickness_plateFins) annotation (Placement(transformation(extent={{-20,-20},{20,20}})));

  ThermofluidStream.Interfaces.Inlet inlet_Cross(redeclare package Medium = Medium_CrossChannel) annotation (
      Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={-100,0}), iconTransformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={-100,0})));
  ThermofluidStream.Interfaces.Outlet outlet_Cross(redeclare package Medium = Medium_CrossChannel) annotation (
      Placement(transformation(extent={{94,-10},{114,10}}), iconTransformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={104,0})));
  ThermofluidStream.Interfaces.Inlet inlet_Pipe(redeclare package Medium = Medium_PipeChannel) annotation (
      Placement(transformation(
        extent={{10,-10},{-10,10}},
        rotation=-90,
        origin={0,-102})));
  ThermofluidStream.Interfaces.Outlet outlet_Pipe(redeclare package Medium = Medium_PipeChannel) annotation (
      Placement(transformation(
        extent={{10,-10},{-10,10}},
        rotation=-90,
        origin={0,100}), iconTransformation(
        extent={{10,-10},{-10,10}},
        rotation=-90,
        origin={0,102})));

  Modelica.Units.SI.MassFlowRate m_flow_cond_total = sum(hxCellPair.m_flow_pipe_cond);

  Modelica.Units.SI.HeatFlowRate Q_flow_A = sum(hxCellPair.Q_flow_cross);
  Modelica.Units.SI.HeatFlowRate Q_flow_B = sum(hxCellPair.Q_flow_pipe);
  Modelica.Units.SI.Mass         Mass_HX =  sum(hxCellPair.TinMass);

  ThermofluidStream.HeatExchangers.Internal.DiscretizedHEXSummary summary "Summary record of Quantities";

  ThermofluidStream.Processes.FlowResistance flowResistanceA[nCells](
    redeclare package Medium = Medium_CrossChannel,
    each r(each displayUnit="mm") = 0.025,
    each l=1,
    redeclare function pLoss = ThermofluidStream.Processes.Internal.FlowResistance.linearQuadraticPressureLoss (
          each k=50)) annotation (Placement(transformation(extent={{30,-10},{50,10}})));
  ThermofluidStream.Topology.JunctionN junctionN(redeclare package Medium = Medium_CrossChannel, N=nCells)
    annotation (Placement(transformation(extent={{60,-10},{80,10}})));
  ThermofluidStream.Topology.SplitterN splitterN(redeclare package Medium = Medium_CrossChannel, N=nCells)
    annotation (Placement(transformation(extent={{-60,-10},{-40,10}})));
protected
  outer ThermofluidStream.DropOfCommons dropOfCommons;

/*  function efficency = ThermofluidStream.HeatExchangers.Internal.calculateEfficency (
  redeclare package Medium_CrossChannel = Medium_CrossChannel,
  redeclare package Medium_PipeChannel = Medium_PipeChannel);*/

initial equation

  if initializeMassFlow then
    inlet_Pipe.m_flow = m_flow_0_B;
    flowResistanceA.m_flow = m_flow_0_A/nCells*ones(nCells);
  else
//    flowResistanceA[1].m_flow = 0;
    for i in 1:nCells-1 loop
      flowResistanceA[i + 1].m_flow = flowResistanceA[1].m_flow;
    end for;
  end if;

//--  hxCellPair[1].conductionElement_pipe.v_mean = 0;

equation
  assert(
    inlet_Cross.m_flow > m_flow_assert,
    "Negative massflow at Air inlet",
    dropOfCommons.assertionLevel);
  assert(
    inlet_Pipe.m_flow > m_flow_assert,
    "Negative massflow at Refigerant inlet",
    dropOfCommons.assertionLevel);

  //Summary record
  summary.Tin_B =Medium_PipeChannel.temperature(inlet_Pipe.state);
  summary.Tin_A =Medium_CrossChannel.temperature(inlet_Cross.state);
  summary.Tout_B =Medium_PipeChannel.temperature(outlet_Pipe.state);
  summary.Tout_A =Medium_CrossChannel.temperature(outlet_Cross.state);
  summary.hin_B =Medium_PipeChannel.specificEnthalpy(inlet_Pipe.state);
  summary.hin_A =Medium_CrossChannel.specificEnthalpy(inlet_Cross.state);
  summary.hout_B =Medium_PipeChannel.specificEnthalpy(outlet_Pipe.state);
  summary.hout_A =Medium_CrossChannel.specificEnthalpy(outlet_Cross.state);
  summary.dT_A = summary.Tout_A - summary.Tin_A;
  summary.dT_B = summary.Tout_B - summary.Tin_B;
  summary.dh_A = summary.hout_A - summary.hin_A;
  summary.dh_B = summary.hout_B - summary.hin_B;
  summary.Q_flow_A=Q_flow_A;
  summary.Q_flow_B=Q_flow_B;
  summary.efficiency = 0;
/*  summary.efficency =if calculate_efficency then efficency(
    inlet_Cross.state,
    inlet_Pipe.state,
    outlet_Cross.state,
    outlet_Pipe.state,
    inlet_Cross.m_flow,
    inlet_Pipe.m_flow,
    Q_flow_A) else 0;*/

  //Connecting equations (to interconnect pipes)

  //Fluid Side B
  for i in 1:nCells-1 loop
  end for;

    for i in 1:nCells-1 loop
    connect(hxCellPair[i].outlet_Pipe, hxCellPair[i + 1].inlet_Pipe);
  end for;

  connect(inlet_Cross, splitterN.inlet)
    annotation (Line(
      points={{-100,0},{-60,0}},
      color={28,108,200},
      thickness=0.5));
  connect(splitterN.outlets, hxCellPair.inlet_Cross)
    annotation (Line(
      points={{-40,0},{-20,0}},
      color={28,108,200},
      thickness=0.5));

  connect(hxCellPair.outlet_Cross, flowResistanceA.inlet)
    annotation (Line(
      points={{20,0},{30,0}},
      color={28,108,200},
      thickness=0.5));
connect(flowResistanceA.outlet, junctionN.inlets) annotation (Line(
      points={{50,0},{60,0}},
      color={28,108,200},
      thickness=0.5));
  connect(junctionN.outlet, outlet_Cross)
    annotation (Line(
      points={{80,0},{104,0}},
      color={28,108,200},
      thickness=0.5));
  connect(inlet_Pipe, hxCellPair[1].inlet_Pipe)
    annotation (Line(
      points={{0,-102},{0,-20}},
      color={28,108,200},
      thickness=0.5));
  connect(hxCellPair[nCells].outlet_Pipe, outlet_Pipe)
    annotation (Line(
      points={{0,20},{0,100}},
      color={28,108,200},
      thickness=0.5));
  annotation (
    Icon(
      graphics={
        Text(
          extent={{-10,76},{2,64}},
          lineColor={28,108,200},
          pattern=LinePattern.Dash,
          textString="..."), Text(
          extent={{22,-78},{80,-92}},
          textColor={0,140,72},
          textString=DynamicSelect("Mass HX",
              String(
              Mass_HX,
              minimumLength=5,
              significantDigits=4) + " kg"))}),
    Documentation(info="<html>
    A general description of cross flow heat exchagers is given in the package.   
    <p>
    <strong>Pipe Flow Media</strong>
    <p>
    1 substance in 1 or 2 phases like refrigerents, water, glycol, thermal oil<br>
    <p>
    <strong>Cross Flow Media</strong>
    <p>
    substance in 1 phase like liquid water, glycol, thermal oil, non-condensing gas<br>
    <p>
</html>"),
    Diagram(graphics={        Text(
          extent={{-58,94},{-18,72}},
          textColor={28,108,200},
          textString="Pipe Flow"), Text(
          extent={{56,-10},{96,-32}},
          textColor={28,108,200},
          textString="Cross Flow")}));
end HX_1PhaseCross_2PhasePipe;
