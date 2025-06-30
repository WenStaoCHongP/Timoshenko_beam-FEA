function Cv=VALUE(W,dispmax,stressmax,frequency,penal)
%----------------------------------------------------------------------------------------
%  compute object value by penalty method
%---------------------------------------------------------------------------------------
if   dispmax<=1e-3
      d1=0.0;
else
      d1=(dispmax-1e-3)/(1e-3); 
end
if   stressmax<=1.732e8
      d2=0.0;
else
      d2=(stressmax-1.732e8)/(1.732e8); 
end
if   frequency>=80
      d3=0.0;
else
      d3=(80-frequency)/80; 
end
   d=d1+d2+d3;
   Cv=1/(W*(1+penal*d));
