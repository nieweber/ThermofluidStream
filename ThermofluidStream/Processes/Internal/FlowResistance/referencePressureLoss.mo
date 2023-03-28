within ThermofluidStream.Processes.Internal.FlowResistance;
function referencePressureLoss "Pressure loss function based on reference values"
  extends Internal.FlowResistance.partialPressureLoss;

  input SI.Pressure dp_ref "Reference pressure drop"
    annotation(Dialog(enable=true));
  input SI.MassFlowRate m_flow_ref "Reference mass flow rate" annotation(Dialog(enable=true));
  input SI.Density rho_ref "Reference density" annotation(Dialog(enable=true));


  input ThermofluidStream.Processes.Internal.ReferencePressureDropFunction dp_function=ThermofluidStream.Processes.Internal.ReferencePressureDropFunction.linear "Pressure drop function"
  annotation (Dialog(enable=true),
  choices(
      choice=ThermofluidStream.Processes.Internal.ReferencePressureDropFunction.linear "Linear",
      choice=ThermofluidStream.Processes.Internal.ReferencePressureDropFunction.quadratic "Quadratic",
      choice=ThermofluidStream.Processes.Internal.ReferencePressureDropFunction.customExponent "Custom exponent"));

  input Real m( unit = "1") = 1.5 "Exponent for pressure drop function"
  annotation(Dialog(enable=(dp_corr == ThermofluidStream.Processes.Internal.ReferencePressureDropFunction.customExponent)));
algorithm

  if dp_function == ThermofluidStream.Processes.Internal.ReferencePressureDropFunction.linear then
    result.dp := rho_ref/rho*dp_ref*m_flow/m_flow_ref;
  elseif dp_function == ThermofluidStream.Processes.Internal.ReferencePressureDropFunction.quadratic then
    result.dp := rho_ref/rho*dp_ref*Modelica.Fluid.Utilities.regSquare(m_flow/m_flow_ref);
  elseif dp_function == ThermofluidStream.Processes.Internal.ReferencePressureDropFunction.customExponent then
    result.dp := rho_ref/rho*dp_ref*Modelica.Fluid.Utilities.regPow(m_flow/m_flow_ref, m);
  end if;

  annotation (Documentation(info="<html>
<p>This function calculates the pressure drop according a reference pressure drop at a reference mass flow rate.</p>
<p>The correlation can be assumed to be linear, quadratic or through a custom exponent.</p>
</html>"));
end referencePressureLoss;
