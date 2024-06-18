within ThermofluidStream.HeatExchangersPhysical.Internal;
function transitionFunction "Gradual Transition between laminar and turbulent phase"
  extends Modelica.Icons.Function;

  input Real y_low  "Value of y at x_low";
  input Real y_high "Value of y at x_high";
  input Real x_low;
  input Real x_high;
  input Real x;
  output Real y;

protected
  Real part_high;

algorithm

  if     noEvent(x < x_low) then y := y_low;
  elseif noEvent(x > x_high) then y := y_high;
  else part_high := (x-x_low)/(x_high-x_low);
       y := y_low*(1-part_high) + y_high*part_high;
  end if;

end transitionFunction;
