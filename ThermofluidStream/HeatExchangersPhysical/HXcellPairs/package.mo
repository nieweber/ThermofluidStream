within ThermofluidStream.HeatExchangersPhysical;
package HXcellPairs "Heat Exchange in Cells between Cross Element and Pipe Element"
  annotation (
    Documentation(info="<html>
    The components in this package model the temperature <b>T<sub>Wall</sub></b> of the wall and fins between one conduction element of the cross channel and one conduction element of the pipe channel.<br>
    The temperature results from the heat flows from/to the fluid in the pipe and to/from the fluid of the cross flow surrounding the pipe.<br>
    <p>
    According to the material (copper, aluminium or steel) the heat capacity of the wall is considered.<br>
    Since the thermal conduction of the wall material is much higher than the thermal conduction of the fluids, the temperature within a cell is assumed to be equal all over the wall and fin area.
    
    <p>
    <img src=\"modelica://ThermoFluidStreamPlus/HeatExchangers/HeatExchangersCross/Resources/HX_HeatFlow_Cell.png\">
    <img src=\"modelica://ThermoFluidStreamPlus/HeatExchangers/HeatExchangersCross/Resources/HX_CellPair.png\">
    <p>

    <p>
    <strong>Cross Flow Media</strong>
    <p>
    substance in 1 phase like liquid water, glycol, thermal oil, non-condensing gas<br>
    Moist Air: condensing or non-condensing<br>
    <p>
    Condensing Air leaves the cross flow channel in 2 flows: moist air and liquid water

    <p>
    <strong>Pipe Flow Media</strong>
    <p>
    1 substance in 1 or 2 phases like refrigerents, water, glycol, thermal oil<br>
    Moist Air: condensing or non-condensing<br>
    <p>
    Condensing Air leaves the pipe in 2 flows: moist air and liquid water

    </html>"));
end HXcellPairs;
