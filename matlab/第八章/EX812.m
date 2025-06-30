% -----------------------------------------------------%   Ex 8.1.2
%   This program is used to optimize cross-section areas of 
%   the space 72-bar truss based on genetic algorithm
% --------------------------------------------------------------------------------
% Initialization of parameters 
% --------------------------------------------------------------------------------
   LENGTH=4;                          %the chromosome length of design variable
   gan=72;                            %the number of gan
   CHROMLENGTH=LENGTH*gan';      %totle length of chromosome
   PopSize=400;                       %population size
   Pc=0.8;                            %probability of crossover
   Pm=0.02;                           %probability of mutation
   penal=100;                        %the initial penalty
 %------------------------------------------------
 % Definition of data structure
 %------------------------------------------------
  for loopi=1:PopSize
     for loopj=1:CHROMLENGTH 
        populationchrom(loopi,loopj)=0.0;
     end
  end
  for loopi=1:3
    for loopj=1:CHROMLENGTH 
        bestworstchrom(loopi,loopj)=0.0;
end 
end
%----------------------------------------------
%  Begin the main program
% ----------------------------------------------
 populationchrom=GenerateInitialPopulation(PopSize,CHROMLENGTH,LENGTH,gan); populationfitness=CalculateObjectValue(populationchrom,LENGTH,PopSize,gan,penal);  bestworstchrom=FindBestAndWorstIndividual(populationchrom,populationfitness,PopSize,bestworstchrom); 
bestworstchrom(3,:)=bestworstchrom(1,:);
%-------the ceasing condition--------------
 itercri=1;
  generation=0;
   while itercri
     populationchrom=GenerateNextPopulation(PopSize,populationchrom,CHROMLENGTH,LENGTH,Pc,Pm,populationfitness);
 populationfitness=CalculateObjectValue(populationchrom,LENGTH,PopSize,gan,penal); bestworstchrom=FindBestAndWorstIndividual(populationchrom,populationfitness,PopSize,bestworstchrom); 
%----output the result of current population----------
     for i=1:3
       area=Decode(bestworstchrom(i,:),LENGTH,gan);
        [W,dispmax,stressmax,frequency]=Trusssizeopt(area);
        C=VALUE(W,dispmax,stressmax,frequency,penal);
        bestworstfitness(i)=C;
    end

   if bestworstfitness(1) > bestworstfitness(3)
        bestworstchrom(3,:)=bestworstchrom(1,:);
   end
     area=Decode(bestworstchrom(3,:),LENGTH,gan);
     [W,dispmax,stressmax,frequency]=Trusssizeopt(area)
     sumfitness=0.0;
     for i=1:PopSize
        sumfitness=sumfitness+populationfitness(i);
     end
     meanfitness=(sumfitness-min(populationfitness)+bestworstfitness(3))/PopSize;
     generation=generation+1
     maxfitdata(generation)=bestworstfitness(3);
     meanfitdata(generation)=meanfitness;
     Edata(generation)=W;
    if generation>50
         h1=abs(Edata(generation)-Edata(generation-50))/Edata(generation);
         if h1<1e-4& frequency>=80&stressmax<=173.2e6&dispmax<=1e-3
             itercri=0;
         else
             itercri=1;
         end
    end
    if penal>=2;
     penal=penal/2;
    end
  end
%------------------------------------------------------------------------------------
% The end
%------------------------------------------------------------------------------------
