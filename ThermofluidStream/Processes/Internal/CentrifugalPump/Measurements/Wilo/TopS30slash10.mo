within ThermofluidStream.Processes.Internal.CentrifugalPump.Measurements.Wilo;
record TopS30slash10 "Pump data for a staged Wilo-Top-S 30/10 pump"
  extends ThermofluidStream.Processes.Internal.CentrifugalPump.Measurements.GenericPump(
    rho_ref=1000,
    w_ref=Modelica.Units.Conversions.from_rpm(2650),
    V_flow={3.51617440225e-06,0.000523909985935,0.000847398030942,0.00115682137834,0.00148382559775,0.00177918424754,0.00206047819972,0.00229254571027,0.0025,0.00271097046414,0.00282700421941},
    P={192.787993617,238.762280675,272.937843988,301.381574494,326.18493197,344.613937245,359.220335761,369.040720135,376.979332273,382.887700535,382.887700535},
    dp={110125.414283,106765.16619,100392.787862,92220.6199738,81651.34693,70879.932776,60307.3708281,50930.4520427,41152.0339559,29575.0912725,22388.3296727},
    final setHead=false,
    final head={0});  // head ist not used since setHead=false. I set final head = {0}; to hide it in the parameter window (for less confusion i hope)
  annotation (
defaultComponentPrefixes="parameter",
defaultComponentName="per",
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
