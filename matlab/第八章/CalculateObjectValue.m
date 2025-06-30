function populationfitness=CalculateObjectValue(populationchrom,LENGTH,PopSize,gan,penal)
%-----------------------------------------------------------------------------------
% to calculate population fitness
%------------------------------------------------  

for loopi=1:PopSize
        area=Decode(populationchrom,LENGTH,gan);
        [W,dispmax,stressmax,frequency]=Trusssizeopt(area);
        Cv=VALUE(W,dispmax,stressmax,frequency,penal); 
        populationfitness(loopi)=Cv;
end
