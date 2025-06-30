function [U5]=reddistrib(nelx,nely,x,dU,dU1,lamda5,lamda6,p,a,b,h)
%-------------------------------------------------------------------------------------
% MESH-INDEPENDENCY FILTER
%-------------------------------------------------------------------------------------
for ely = 1:nely
    for elx = 1:nelx
      U4(ely,elx)=((lamda5/(4*p*a*b))*dU(ely,elx)*x(ely,elx)+(lamda6/(4*p*a*b))*dU1(ely,elx)*x(ely,elx))^h;
  end
end
for ely = 1:nely+1
    for elx = 1:nelx+1
        sum1=0;
        MM1=0;
        aa1=[elx ely;elx-1 ely-1;elx-1 ely;elx ely-1];
        for i=1:4
         if(aa1(i,:)>0&aa1(i,1)<=nelx&aa1(i,2)<=nely)
          sum1=sum1+U4(aa1(i,2),aa1(i,1));
           MM1=MM1+1;
         end
        end
       U6(ely,elx)=sum1/MM1;
    end
end
for ely = 1:nely
    for elx = 1:nelx
     sum11=0;
     MM11=0;
     aa11=[elx ely;elx+1 ely+1;elx ely+1;elx+1 ely];
      for i=1:4
       if(aa11(i,:)>0&aa11(i,1)<=nelx+1&aa11(i,2)<=nely+1)
       sum11=sum11+U6(aa11(i,2),aa11(i,1));
       MM11=MM11+1;
       end
      end
     U5(ely,elx)=sum11/MM11;
    end
end
