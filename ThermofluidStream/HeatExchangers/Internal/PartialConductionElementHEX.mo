within ThermofluidStream.HeatExchangers.Internal;
partial model PartialConductionElementHEX "Parent for CEs for discretizedHEX"
  extends Processes.Internal.PartialConductionElement(
    final init= Processes.Internal.InitializationMethodsCondElement.inlet,
    final neglectPressureChanges=true);

    parameter SI.Area A = 1 "Contact area of volume with medium";

    parameter Integer nCellsParallel = 1 "Number of parallel discretization elements";

    constant SI.CoefficientOfHeatTransfer U_min = 1 "Minimum heat transfer coefficient for temperature adaption at zero massflow";

    SI.CoefficientOfHeatTransfer U "Heat transfer coefficient to medium";

    Modelica.Blocks.Interfaces.RealInput DT1_input(unit = "K")
      annotation (Placement(transformation(extent={{-20,-20},{20,20}},
        rotation=90,
        origin={-50,-100})));
    Modelica.Blocks.Interfaces.RealInput DT2_input(unit = "K")
      annotation (Placement(transformation(extent={{-20,-20},{20,20}},
        rotation=90,
        origin={50,-100})));

equation
  k = U*A;

  DT1 = DT1_input;
  DT2 = DT2_input;

  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
end PartialConductionElementHEX;
