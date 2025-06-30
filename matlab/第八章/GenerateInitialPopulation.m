function populationchrom=GenerateInitialPopulation(PopSize,CHROMLENGTH,LENGTH,gan)
%-----------------------------------------------------------------------------
% generate the first population
%-----------------------------------------------------------------------------

  for loopi=1:PopSize
      for loopj=1:(LENGTH*gan)
        if rand <=0.5
           populationchrom(loopi,loopj)=1;
        else
           populationchrom(loopi,loopj)=0;
        end
     end
  end
