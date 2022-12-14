within ThermofluidStream.Utilities.Functions;
function logMeanRegularized

  import Modelica.Math;

  input SI.Temperature DT1;
  input SI.Temperature DT2;
  //input Real gamma = 15;
  input Boolean applySmoothing;
  output SI.Temperature DTLM;

protected
  parameter Integer N = 4 "Order of quadrature";

  Real f;

  Real grid[N+1] = {0, 1/8, 1/4, 1/2, 1};
  Real dx[N] = zeros(N);


  Real DT1a;
  Real DT2a;

  Real DT1s "DT1 for smoothing";
  Real DT1_norm = 1e-5 "Normalization for both inputs being zero";

  parameter Real alpha_gradient = 0.0;
  parameter Real gamma = 15;


algorithm

   if abs(DT1) > abs(DT2) then
     DT1a := DT1;
     DT2a := DT2;
   else
     DT1a := DT2;
     DT2a := DT1;
     end if;

   DTLM := 0;

   f := abs(DT2a*DT1a/(DT1a^2+DT1_norm^2));

//    //Loop to apply trapezoidal rule
//    for i in 2:N+1 loop
//
//      dx[i-1] := grid[i] - grid[i-1];
//
//      DTLM := DTLM + (1/2)*(f^(grid[i-1])+f^(grid[i]))*dx[i-1];
//
//    end for;

    DTLM := (1.0/4)*(f^(1.0/2)+f);
    DTLM := DTLM + (1.0/8)*(f^(1.0/4)+f^(1.0/2));
    DTLM := DTLM + (1.0/16)*(f^(1.0/8)+f^(1.0/4));
    DTLM := DTLM + (1.0/16)*(f^(1.0/8)+1);

   //Smoothing the edges
   if DT2a > 0 then

     DT1s :=min(DT1a, gamma*DT2a);
   else
     DT1s :=max(DT1a, gamma*DT2a);

      end if;

   if applySmoothing then
     DTLM := DT1s*DTLM;
   else
     DTLM := DT1a*DTLM;
   end if;

   if DT1*DT2 <=0 then

     DTLM := 0;
     return;
   end if;

   if abs(DT2) - abs(DT1) == 0 then

     DTLM := (DT1+DT2)/2;
   end if;

   //"Tilting" the planes in the area where x*L<=0
   DTLM := (1 - alpha_gradient)*DTLM + alpha_gradient*(DT1+DT2)/2;

end logMeanRegularized;
