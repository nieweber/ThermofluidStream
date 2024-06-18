within ThermofluidStream.HeatExchangersPhysical.Internal;
function portion "Gradual Portion between laminar and turbulent phase"
  extends Modelica.Icons.Function;

  input  Real x_low;
  input  Real x_high;
  input  Real x;
  output Real y_low;
  output Real y_high;

algorithm

  y_high := if noEvent(x < x_low)  then 0
        elseif noEvent(x > x_high) then 1
                                   else (x-x_low)/(x_high-x_low);

  y_low := 1 - y_high;

end portion;
