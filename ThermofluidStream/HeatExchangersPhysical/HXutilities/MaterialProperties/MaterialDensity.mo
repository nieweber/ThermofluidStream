within ThermofluidStream.HeatExchangersPhysical.HXutilities.MaterialProperties;
function MaterialDensity "Density of a solid material"

  import Modelica.Units.SI;

  input Material material;
  output SI.Density density;

  constant SI.Density  MaterialDensity[3] = {2710, 8940, 7880};

// https://de.wikibooks.org/wiki/Tabellensammlung_Chemie/_Dichte_fester_Stoffe

algorithm

  density :=MaterialDensity[material];

end MaterialDensity;
