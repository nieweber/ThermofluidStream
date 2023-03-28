within ThermofluidStream.Processes.Internal.FlowResistance;
partial function partialPressureLoss
  "Interface for pressure loss functions"
  input SI.MassFlowRate m_flow "Mass flow rate";
  input SI.Density rho "Medium density at inlet";
  input SI.Density rho_out "Medium density at outlet";
  input SI.DynamicViscosity mu "Medium dynamic viscosity";
  input SI.Pressure p_out "Pressure at outlet";
  input SI.Temperature T "Temperature";
  input SI.SpecificHeatCapacity Rs "Specific gas constant";

  output dpRecord result "Collects results";

  //output SI.Pressure pressureLoss "Pressure lost in Pipe";
  //output SI.Area crossSectionArea "Effective cross section area";
  //SI.Length d_h "Hydraulic diameter of resistance";
  //output SI.Velocity v_mean "Mean flow velocity";
protected
  SI.Length r_h "Hydraulic radius";

  SI.Length d_h "Hydraulic radius";

  annotation(Inline=true, smoothOrder=100,
    Documentation(info="<html>
<p>Interface definition for a pressure loss in a pipe. Inputs are information about flow condition and the medium as well as the geometry of the pipe, output is the pressure drop.</p>
</html>"));
end partialPressureLoss;
