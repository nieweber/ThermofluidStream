within ThermofluidStream.Processes.Internal.FlowResistance;
function zetaPressureLoss "Pressure loss function based on zeta value"
  extends Internal.FlowResistance.partialPressureLoss;

  input Real zeta( unit = "1") "Zeta value of component"
    annotation(Dialog(enable=true));
  input Boolean areaFromRadius = true "Calculate reference area from pipe radius?"
    annotation(Dialog(enable = true));
  input SI.Area A = Modelica.Constants.pi*r*r "Reference area from parameter"
    annotation(Dialog(enable=(not areaFromRadius)));

protected
  SI.Area A_zeta "Reference area either from radius or set by parameter";

algorithm

  if areaFromRadius then
    A_zeta :=Modelica.Constants.pi*r*r;
  else
    A_zeta :=A;
  end if;

  pressureLoss :=zeta/(2*rho)*Modelica.Fluid.Utilities.regSquare(m_flow/A_zeta);

  annotation (Documentation(info="<html>
<p>For specific components (armatures, fittings, pipe sections, grids, ...), the zeta value is often given in the data sheet.</p>
<p>Together with a given reference area A, the pressure drop can be calculated:</p>
<p style=\"margin-left: 40px;\"><span style=\"font-family: Courier New;\">dp := zeta/(2*rho)*m_flow^2/A^2</span></p>
<p>The square of the mass-flow is regularized using the regSquare function from the MSL.</p>
</html>"));
end zetaPressureLoss;
