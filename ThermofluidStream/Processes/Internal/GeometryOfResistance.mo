within ThermofluidStream.Processes.Internal;
type GeometryOfResistance = enumeration(
    circular,
    rectangle,
    other) "Provides the choice of different geometries"                  annotation (
    Documentation(info="<html>
<p>Enumeration of materials where the roughness is approximately known. Used for a dropdown of selectable materials.</p>
</html>"));
