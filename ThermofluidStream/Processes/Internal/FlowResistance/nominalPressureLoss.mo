within ThermofluidStream.Processes.Internal.FlowResistance;
function nominalPressureLoss
  "Pressure loss function based on nominal values"
  extends Internal.FlowResistance.partialPressureLoss;

  input SI.Pressure dp_ref "Nominal pressure drop"
    annotation(Dialog(enable=true));
  input SI.MassFlowRate m_flow_ref "Nominal mass flow rate" annotation(Dialog(enable=true));

  input ThermofluidStream.Processes.Internal.PressureDropCorrelation dp_corr = ThermofluidStream.Processes.Internal.PressureDropCorrelation.linear
  "Pressure drop correlation"
    annotation (Dialog(enable=true),
     choices(
      choice=ThermofluidStream.Processes.Internal.PressureDropCorrelation.linear "Linear",
      choice=ThermofluidStream.Processes.Internal.PressureDropCorrelation.quadratic "Quadratic",
      choice=ThermofluidStream.Processes.Internal.PressureDropCorrelation.customExponent "Custom exponent"));

  input Real m( unit = "1") = 1.5 "Exponent for pressure drop correlation"
  annotation(Dialog(enable = (dp_corr == ThermofluidStream.Processes.Internal.PressureDropCorrelation.customExponent)));
algorithm

  if dp_corr == ThermofluidStream.Processes.Internal.PressureDropCorrelation.linear then
    pressureLoss :=dp_ref*m_flow/m_flow_ref;
  elseif dp_corr == ThermofluidStream.Processes.Internal.PressureDropCorrelation.quadratic then
    pressureLoss := dp_ref*Modelica.Fluid.Utilities.regSquare(m_flow/m_flow_ref);
  elseif dp_corr == ThermofluidStream.Processes.Internal.PressureDropCorrelation.customExponent then
    pressureLoss :=dp_ref*Modelica.Fluid.Utilities.regPow(m_flow/m_flow_ref, m);
  end if;

  annotation (Documentation(info="<html>
<p>This function calculates the pressure drop according a reference pressure drop at a reference mass flow rate.</p>
<p>The correlation can be assumed to be linear, quadratic or through a custom exponent.</p>
</html>"));
end nominalPressureLoss;
