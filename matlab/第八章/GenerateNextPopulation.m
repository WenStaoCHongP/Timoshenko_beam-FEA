function populationchrom=GenerateNextPopulation(PopSize,populationchrom,CHROMLENGTH,LENGTH,Pc,Pm,populationfitness)
%-------------------------------------------------------------------------------------
%   generate next population:select,crossover and mutate
%-------------------------------------------------------
%  to reproduce a chromosome by proportional selection
%-------------------------------------------------------
  p=0.0;
  sum=0.0;
  %------------------------------
  %   calculate relative fitness
  %------------------------------
  for i=1:PopSize
     sum=sum+populationfitness(i);
  end
  for i=1:PopSize
     cfitness(i)=populationfitness(i)/sum;
  end
  %-----------------------------
  % calculate cumulative fitness
  %-----------------------------
  for i=2:PopSize
     cfitness(i)=cfitness(i-1)+cfitness(i);
  end
  %----------------------------
  % selection operation
  %----------------------------
  for i=1:PopSize
     p=rand;
     index=1;
     while p > cfitness(index)
        index=index+1;
     end
     newpopulationchrom(i,:)=populationchrom(index,:);
  end
  for i=1:PopSize
     populationchrom(i,:)=newpopulationchrom(i,:);
  end
  
  %--------------------------------------------------------------
  % crossover two chromosome by means of two-point crossover
  %--------------------------------------------------------------
  % make a pair of individual randomly
  %------------------------------------
  for i=1:PopSize
     index(i)=i;
  end
  for i=0:PopSize-1
     point=round(random('unif',0,PopSize-i));
     temp=index(i+1);
     if (point+i+1) > PopSize
        index(i+1)=index(point+i);
        index(point+i)=temp;
     else
        index(i+1)=index(point+i+1);
        index(point+i+1)=temp;
     end
  end
  %----------------------------------
  % two point crossover operation
  %---------------------------------
  for i=1:2:PopSize
     p=rand;
     if p<Pc
        point(1)=round(random('unif',1,CHROMLENGTH));
        point(2)=round(random('unif',1,CHROMLENGTH));
        if point(1)==point(2)
          ch=populationchrom(index(i),1:point(1));
          populationchrom(index(i),1:point(1))=populationchrom(index(i+1),1:point(1));
          populationchrom(index(i+1),1:point(1))=ch;
        elseif point(1)>point(2)
          ch=populationchrom(index(i),point(2):point(1));
          populationchrom(index(i),point(2):point(1))=populationchrom(index(i+1),point(2):point(1));
          populationchrom(index(i+1),point(2):point(1))=ch;
        else
          ch=populationchrom(index(i),point(1):point(2));
          populationchrom(index(i),point(1):point(2))=populationchrom(index(i+1),point(1):point(2));
          populationchrom(index(i+1),point(1):point(2))=ch;
       end
     end
     end
  %--------------------------------------------------------------
  %           mutation of a chromosome
  %-------------------------------------------------------------
  for i=1:PopSize
     for j=1:CHROMLENGTH
        p=rand;
        if p<Pm
           if populationchrom(i,j)==0;
              populationchrom(i,j)=1;
           else
              populationchrom(i,j)=0;
           end
        end
     end
  end
