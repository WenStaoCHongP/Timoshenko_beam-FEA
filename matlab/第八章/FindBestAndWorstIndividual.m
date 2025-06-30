function bestworstchrom=FindBestAndWorstIndividual(populationchrom,populationfitness,PopSize,bestworstchrom); 
%-------------------------------------------------------------------------------------------------
% to find out the best and worst individual so far current generation
%-------------------------------------------------------------------------------------
  bestworstfitness(1)=populationfitness(1);
  bestworstfitness(2)=populationfitness(1);
   flag=1;
 while flag
       num=1;
    for loopi=2:PopSize
        num=num+1;
     if  populationfitness(loopi)>=bestworstfitness(1)
        bestworstchrom(1,:)=populationchrom(loopi,:);
        bestworstfitness(1)=populationfitness(loopi);
     end
     if  populationfitness(loopi)<=bestworstfitness(2)
        bestworstchrom(2,:)=populationchrom(loopi,:);
        bestworstfitness(2)=populationfitness(loopi);
     end
     if num==PopSize
         flag=0;
     end
  end
 end
