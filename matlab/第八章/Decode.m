function area=Decode(populationchrom,LENGTH,gan)
%----------------------------------------------------------------------------------
%  decode
%  ouput variables:
%   area =cross sectional areas of structural elements
%----------------------------------------------------------------------------------
% extract structural cross section areas
%-----------------------------------------
  vdarea=1e-6*[3125  80 100 125 180   259.2  320 405 500  672.2  793.8  1125  1280  1620  2000  2420];
  area=zeros(gan(1),1);
  for loopi=1:gan(1)
      temp=0;
     for loopj=1:LENGTH(1)
        if populationchrom((loopi-1)*LENGTH(1)+loopj)==1
           temp=temp+2^(loopj-1);
        end
     end
     temp=temp+1;
     area(loopi)=vdarea(temp);
  end
