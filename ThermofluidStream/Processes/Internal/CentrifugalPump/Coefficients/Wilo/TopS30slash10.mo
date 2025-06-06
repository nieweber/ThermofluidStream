within ThermofluidStream.Processes.Internal.CentrifugalPump.Coefficients.Wilo;
record TopS30slash10 "Pump coefficients of a staged Wilo-Top-S 30/10 pump"
  extends ThermofluidStream.Processes.Internal.CentrifugalPump.Coefficients.GenericPump
                                                                               (
    final setHead=true,
    final dp_peak=0,
    w_ref=Modelica.Units.Conversions.from_rpm(2650),
    rho_ref=1000,
    V_flow_peak=0.001591,
    head_peak=8.04,
    eta_peak=0.378,
    V_0n=2.018,
    h_0n=1.409,
    P_0n=0.568);

  annotation (
defaultComponentPrefixes="parameter",
Documentation(info="<html>
<p>
<a href=\"modelica://ThermofluidStream/Resources/Doku/Pumps/WiloTopS30_10_data.pdf\">Data sheet of Wilo-Top-S 30/10</a>
</p>
<p>
<img src=\"modelica://ThermofluidStream/Resources/Doku/Pumps/WiloTopS30_10_head.png\" alt=\"Wilo-Top-S 30/10 head\" width=\"500\">
<img src=\"modelica://ThermofluidStream/Resources/Doku/Pumps/WiloTopS30_10_power.png\" alt=\"Wilo-Top-S 30/10 power\" width=\"500\">
<img src=\"modelica://ThermofluidStream/Resources/Doku/Pumps/WiloTopS30_10_efficiency.png\" alt=\"Wilo-Top-S 30/10 efficiency\" width=\"500\">
</p>
</html>",
  revisions="<html>
<ul>
<li>
June 17, 2024, by Raphael Gebhart:
<br>
Transferred from <a href=\"https://github.com/lbl-srg/modelica-buildings\">Buildings Library</a> to  <a href=\"https://github.com/DLR-SR/ThermofluidStream\">ThermofluidStream Library</a>.
</li>
<li>
March 29, 2023, by Hongxiang Fu:<br>
Deleted angular speed parameters with the unit rpm.
This is for
<a href=\"https://github.com/ibpsa/modelica-ibpsa/issues/1704\">IBPSA, #1704</a>.
</li>
<li>
October 14, 2021, by Hongxiang Fu:<br>
Rewrote the statements using <code>use_powerCharacteristic</code>
to support the implementation of
<a href=\"Modelica://Buildings.Fluid.Movers.BaseClasses.Euler\">
<code>Buildings.Fluid.Movers.BaseClasses.Euler</code></a>.
This is for
<a href=\"https://github.com/lbl-srg/modelica-buildings/issues/2668\">#2668</a>.
</li>
<li>
February 17, 2016, by Michael Wetter:<br>
Updated parameter names for
<a href=\"https://github.com/ibpsa/modelica-ibpsa/issues/396\">#396</a>.
</li>
<li>
February 17, 2016, by Michael Wetter:<br>
Updated parameter names for
<a href=\"https://github.com/ibpsa/modelica-ibpsa/issues/396\">#396</a>.
</li>
<li>
January 19, 2016, by Filip Jorissen:<br>
First implementation.
</li>
</ul>
</html>"));
end TopS30slash10;
