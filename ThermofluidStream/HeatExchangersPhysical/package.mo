within ThermofluidStream;
package HeatExchangersPhysical "Cross Flow Heat Exchanger Package - Version 1.0.0"

  annotation (
    versionDate="2024-06-13",
    Documentation(info="<html>
<p>
<strong>Motivation</strong>
<p>
This package has been designed to meet the requirements of climate systems, particularly environmental conditioning systems of aircraft.<br>
For this reason the heat exchangers are carried out in cross flow and consider phase changes of refrigerents as well as of humid air.<br>
Condensors, evaporators, air-to-air and air-to-liquid can be modelled.<br>
The behavior of an heat exchanger is derived from geometric data. No coefficients like heat exchange rates are needed.

<p>
<strong>Structure of the Heat Exchangers</strong>
<p>
The heat exchangers work in cross flow. One medium is carried in a meandering pipe, which passes the other medium in a perpendicular flow.<br>
The cross air passage is rectangular with dimensions <b>width_HX</b> and <b>height_HX</b> and is of depth <b>length_HX</b>. <br>
The pipe runs  in horizontal direction through the air passage meander-shaped in <b>nRows</b> rows. <br>
Each row may be divided in <b>nCellPerRow</b> cells and may carry <b>nPlateFinsPerRow</b> fins to enlarge the heat exchange area. <br>
The thickness of the pipe and the plate fins as well as the diameter of the pipe has to be specified.<br>
The heat is exchanged through the pipe walls and the fins. Their material may be copper, aluminium or steel.<br>
Depending on possible flow media 3 different types of heat exchangers are realized.<br>

<p>
<img src=\"modelica://ThermoFluidStreamPlus/HeatExchangers/HeatExchangersCross/Resources/HX_CrossFlow.png\">
<p>

<p>
<strong>Geometric Data</strong>
<p>
Besides the material, the following geometric data define the behavior of a heat exchanger compeletely:
<p>
<table> <tbody>
<tr>   <td>     nRows                                        </td>  <td>   Number of rows                                </td> </tr>
<tr>   <td>     nCellPerRow                                  </td>  <td>   Number of cells per row                       </td> </tr>
<tr>   <td>     nPlateFinsPerRow                             </td>  <td>   Number of plate fins per row                  </td> </tr>
<tr>   <td>     length_HX                                    </td>  <td>   length of heat exchanger = length of fins [m] </td> </tr>
<tr>   <td>     width_HX                                     </td>  <td>   width of heat exchanger [m]                   </td> </tr>
<tr>   <td>     height_HX                                    </td>  <td>   height of heat exchanger [m]                  </td> </tr>
<tr>   <td>     inner_diameter_pipe                          </td>  <td>   inner diameter of pipe channel [m]            </td> </tr>
<tr>   <td>     thickness_pipe                               </td>  <td>   thickness of pipe wall [m]                    </td> </tr>
<tr>   <td>     nCellPerRow                                  </td>  <td>   Number of cells per row                       </td> </tr>
<tr>   <td>     thickness_plateFins &nbsp;&nbsp;&nbsp;&nbsp; </td>  <td>   thickness of plate fins [m]                   </td> </tr>
</tbody> </table>

<p>
<strong>Discretization</strong>
<p>
The geometry is discretized in such a way that each element of a pipe is assigned to an element of the air passage separated by a piece of pipe wall.
The heat transfer takes place between such a pair of conduction elements.</p>
<p>
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

<p>
<strong>Flow Resistance</strong>
<p>
With the given geometric data the flow resistance is calculated, too.<br>

<p>
<strong>Modeling Approach and Validation</strong>
<p>
Modeling was carried out according to the text book from H.D. Baehr and K. Stephan: \"W&auml;rme- und Stoff&uuml;bertragung\".<br>
Condenser and Evaporator cooled with dry air could be validated on a test rig.<br>
The other cases, particulary condensation of humid air could not be tested, because there was no test equipment yet available.

</html>"));
end HeatExchangersPhysical;
