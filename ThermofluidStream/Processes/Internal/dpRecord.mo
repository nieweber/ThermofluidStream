within ThermofluidStream.Processes.Internal;
record dpRecord

  SI.Pressure dp "Pressure drop";
  SI.Length d_h "Hydraulic diameter";
  SI.Area A "Effective area";
  SI.Velocity v_mean "Mean flow velocity";
  Real zeta "Zeta value";

  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
end dpRecord;
