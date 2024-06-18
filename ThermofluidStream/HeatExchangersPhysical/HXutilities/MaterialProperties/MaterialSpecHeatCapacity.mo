within ThermofluidStream.HeatExchangersPhysical.HXutilities.MaterialProperties;
function MaterialSpecHeatCapacity "Specific Heat Capacity of a solid material"

  import Modelica.Units.SI;

  input Material material;
  output SI.SpecificHeatCapacity specHeatCapacity;

  constant SI.SpecificHeatCapacity  MaterialHeatCapacity[3] = { 895, 381, 477};

//  https://de.wikibooks.org/wiki/Tabellensammlung_Chemie/_spezifische_Wärmekapazitäten

algorithm

  specHeatCapacity := MaterialHeatCapacity[material];

end MaterialSpecHeatCapacity;
