function [xnew]=OC(nelx,nely,x,h,U5,q,G,freq,d,g,freq1,p,a,b,dU,G1,G2)
%-------------------------------------------------------------------------------------
% OPTIMALITY CRITERIA UPDATE
%-------------------------------------------------------------------------------------
step1=1;
for ely = 1:nely
   for elx = 1:nelx
     xnew1(ely,elx)=real(x(ely,elx)*U5(ely,elx))^(1/(1+2*h));   xnew2(ely,elx)=real((((freq1/(freq(1)))^d)*(((G1/G2)*(G(ely,elx)/((1/dU(ely,elx)))))^g))*x(ely,elx));
     TT1((elx-1)*nely+ely)=xnew1(ely,elx);
     TT2((elx-1)*nely+ely)=xnew2(ely,elx);
   end
end
 

for ely = 1:nely
    for elx = 1:nelx
     xnew2(ely,elx)=xnew2(ely,elx)*(median(TT1)/median(TT2));
     xnew(ely,elx)=min(max(xnew1(ely,elx),xnew2(ely,elx)),(xnew1(ely,elx)+xnew2(ely,elx))/2);    
     T((elx-1)*nely+ely)=xnew(ely,elx);
   end
end
 Mt=sort(T);
     for i=1:(nely*nelx)
       if((0.00006)<Mt(i)&Mt(i)<(0.006))
        Maa(step1)=Mt(i);
        step1=step1+1;
       end
     end
TT=2*floor(step1*0.01)
for ely = 1:nely
    for elx = 1:nelx
     if((xnew(ely,elx)>0.00006)&(xnew(ely,elx)<0.006)&TT~=0)
        if(xnew(ely,elx)<=Maa(TT))
         xnew(ely,elx)=0.00006;
        end
     end
     TTT((elx-1)*nely+ely)=xnew(ely,elx);   
    end
 end        
for ely = 1:nely
    for elx = 1:nelx
     w(ely,elx)=q*sign(xnew(ely,elx)-mean(TTT))*xnew(ely,elx);
     xnew(ely,elx)=xnew(ely,elx)+w(ely,elx);
     xnew(ely,elx)=min(max(xnew(ely,elx),0.00001),0.006);
     if((abs(x(ely,elx)-0.00001)<=1e-5)|(abs(x(ely,elx)-0.006)<=1e-5))
     xnew(ely,elx)=x(ely,elx);
     end
    end
 end
