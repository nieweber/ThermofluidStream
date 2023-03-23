within ThermofluidStream.Processes.Internal;
type ReferencePressureDropFunction = enumeration(
    linear,
    quadratic,
    customExponent) "Provides the choice of linear or quadratic pressure loss functions"     annotation (
    Documentation(info="<html>
<p>Enumeration of materials where the roughness is approximately known. Used for a dropdown of selectable materials.</p>
</html>"));
