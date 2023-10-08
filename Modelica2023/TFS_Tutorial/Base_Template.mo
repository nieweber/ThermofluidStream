within TFS_Tutorial;
model Base_Template

  extends Modelica.Icons.Example;

  //Parameters
  parameter Modelica.Units.SI.Temperature heatingThresholdLow=283.15   "Lower threshold for active battery heating" annotation(Dialog(group = "Battery"));
  parameter Modelica.Units.SI.Temperature heatingThresholdHigh=298.15       "Upper threshold for active battery heating" annotation(Dialog(group = "Battery"));
  parameter Modelica.Units.SI.Temperature coolingThreshold=303.15         "Threshold for active cooling" annotation(Dialog(group = "Battery"));
  parameter Modelica.Units.SI.Temperature coolingThresholdVCS=313.15         "Threshold for active VCS cooling" annotation(Dialog(group = "Battery"));
  parameter Modelica.Units.SI.Temperature batteryTargetTemperature=298.15 "Desired battery temperature" annotation(Dialog(group = "Battery"));
  parameter Modelica.Units.SI.Temperature T_amb=288.15 "Ambient temperature" annotation(Dialog(group = "Ambient conditions"));

//Media definitions
  package Coolant =
      ThermofluidStream.Media.myMedia.Incompressible.Examples.Glycol47;
  package Air = ThermofluidStream.Media.myMedia.Air.MoistAir;
  package Refrigerant = ThermofluidStream.Media.XRGMedia.R134a_ph;

  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-260,-160},{260,160}})));
end Base_Template;
