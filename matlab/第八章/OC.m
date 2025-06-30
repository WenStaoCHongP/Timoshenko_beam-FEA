function [xnew]=OC(nelx,nely,x,h,U5,q,G,freq,d,g,freq1,p,a,b,dU,G1,G2) 
%--------------------------------------------------------------------------------------
%  OPTIMALITY CRITERIA UPDATE
%-------------------------------------------------------------------------------------
step1=1;
for ely = 1:nely
    for elx = 1:nelx
     xnew1(ely,elx)=real(x(ely,elx)*U5(ely,elx))^(1/(1+2*h));     %  update design variables based displacement sensitivity
     xnew2(ely,elx)=real((((freq1/(freq(1)))^d)*(((G1/G2)*(G(ely,elx)/((1/dU(ely,elx)))))^g))*x(ely,elx));  %  update design variables based frequency sensitivity
     TT1((elx-1)*nely+ely)=xnew1(ely,elx);
     TT2((elx-1)*nely+ely)=xnew2(ely,elx);
   end
end
 

for ely = 1:nely
    for elx = 1:nelx
     xnew2(ely,elx)=xnew2(ely,elx)*(median(TT1)/median(TT2));  %  tension of different update design variables based different sensitivities  
     xnew(ely,elx)=min(max(xnew1(ely,elx),xnew2(ely,elx)),(xnew1(ely,elx)+xnew2(ely,elx))/2);  %  tradeoff of update design variables 
     T((elx-1)*nely+ely)=xnew(ely,elx);
   end
end

%  deletion of certain percentage of elements using Evolutionary Structural Optimization(ESO)
 Mt=sort(T);
     for i=1:(nely*nelx)
       if((0.00005)<Mt(i)&Mt(i)<(0.009))
        Maa(step1)=Mt(i);
        step1=step1+1;
       end
     end
TT=2*floor(step1*0.005)
for ely = 1:nely
    for elx = 1:nelx
     if((xnew(ely,elx)>0.000009)&(xnew(ely,elx)<0.009)&TT~=0)
        if(xnew(ely,elx)<=Maa(TT))
         xnew(ely,elx)=0.000009;
        end
     end
     TTT((elx-1)*nely+ely)=xnew(ely,elx);   
    end
 end        

%  division of design variable based average value
for ely = 1:nely
    for elx = 1:nelx
     w(ely,elx)=q*sign(xnew(ely,elx)-mean(TTT))*xnew(ely,elx);
     xnew(ely,elx)=xnew(ely,elx)+w(ely,elx);
     xnew(ely,elx)=min(max(xnew(ely,elx),0.000009),0.009);
     if((abs(x(ely,elx)-0.000009)<=1e-5)|(abs(x(ely,elx)-0.009)<=1e-5))
     xnew(ely,elx)=x(ely,elx);
     end
    end
 end
