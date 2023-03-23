within ThermofluidStream.Processes.Internal;
type ShapeOfResistance = enumeration(
    circular,
    rectangle,
    other) "Provides the choice of different shapes for cross section"                  annotation (
    Documentation(info="<html>
<p>Enumeration of materials where the roughness is approximately known. Used for a dropdown of selectable materials.</p>
</html>"));
