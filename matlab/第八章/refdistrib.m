function [G,G1,G2]=refdistrib(nelx,nely,G,p,a,b,dU)
%-------------------------------------------------------------------------
% frequency sensitivity and sensitivity redistrubution
%-------------------------------------------------------------------------
G1=0;
G2=0;
for ely = 1:nely+1
    for elx = 1:nelx+1
        sum2=0;
        MM2=0;
        aa2=[elx ely;elx-1 ely-1;elx-1 ely;elx ely-1];
        for i=1:4
         if(aa2(i,:)>0&aa2(i,1)<=nelx&aa2(i,2)<=nely)
          sum2=sum2+G(aa2(i,2),aa2(i,1));
           MM2=MM2+2;
         end
        end
       G3(ely,elx)=sum2/MM2;
    end
end
for ely = 1:nely
    for elx = 1:nelx
     sum21=0;
     MM21=0;
     aa21=[elx ely;elx+1 ely+1;elx ely+1;elx+1 ely];
      for i=1:4
       if(aa21(i,:)>0&aa21(i,1)<=nelx+1&aa21(i,2)<=nely+1)
       sum21=sum21+G3(aa21(i,2),aa21(i,1));
       MM21=MM21+1;
       end
      end
     G(ely,elx)=sum21/MM21;
    end
end
for ely = 1:nely
    for elx = 1:nelx
    G1=G1+(1/dU(ely,elx))*G(ely,elx);
    G2=G2+G(ely,elx)^2;
    end
end
