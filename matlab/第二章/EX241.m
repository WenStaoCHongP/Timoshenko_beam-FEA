%-------------------------------------------------------------------------------------
%  Ex.2.1 Plane 61-bar truss
%-------------------------------------------------------------------------------------
%  To compute the first nine eigenvalues of plane 61-bar truss
% 
E=2.1e11;
A=1e-4;
density=7.3e3; 
node_number=26;
elment_number=61;
nc=zeros(26,2);
nc(1,1)=0;nc(1,2)=0;
nc(2,1)=0;nc(2,2)=1;
for ig=1:12
    for jg=1:2
        nc(ig*2+jg,1)=nc(jg,1)+ig;           %node_coordinate
        nc(ig*2+jg,2)=nc(jg,2);
    end
end
en=zeros(61,2);
en(1,1)=1;    en(1,2)=2;                      %en:element_node  
en(2,1)=1;    en(2,2)=4;
en(3,1)=2;    en(3,2)=3;
en(4,1)=1;    en(4,2)=3;
en(5,1)=2;    en(5,2)=4;
for i=1:11
   for j=1:5
      en(i*5+j,1)=en(j,1)+2*i;
      en(i*5+j,2)=en(j,2)+2*i;
   end
end
en(61,1)=25;    en(61,2)=26;
ed(1:node_number,1:2)=1;                    %ed:elment_displacement 
constraint=[1,1;1,2;25,2];
for loopi=1:length(constraint);
    ed(constraint(loopi,1),constraint(loopi,2))=0;
end
dof=0;
for loopi=1:node_number
     for loopj=1:2
        if ed(loopi,loopj)~=0
            dof=dof+1;
            ed(loopi,loopj)=dof;
        end
    end
end
ek=E*A*[1 0 -1 0;0 0 0 0;-1 0 1 0;0 0 0 0];    %elment_stiffness_matrix
em=(density*A)/6*[2 0 1 0; ...              %elment_mass_matrix  
                   0 2 0 1; ...
                   1 0 2 0; ...
                  0 1 0 2];
k(1:dof,1:dof)=0;                          %structural_stiffness_matrix  
m=k;                                     %structural_mass_matrix,same size with k 
theta(1:61)=0;
el(1:61)=0;
e2s(1:4)=0;                         %  index of transform the elment displament number to structure
for loopi=1:elment_number
    for zi=1:2
        e2s((zi-1)*2+1)=ed(en(loopi,zi),1);
        e2s((zi-1)*2+2)=ed(en(loopi,zi),2);
    end
    el(loopi)=sqrt((nc(en(loopi,1),1)-nc(en(loopi,2),1))^2+(nc(en(loopi,1),2)-nc(en(loopi,2),2))^2);
    theta(loopi)=asin((nc(en(loopi,1),2)-nc(en(loopi,2),2))/el(loopi));
    lmd=[cos(theta(loopi)) sin(theta(loopi)); -sin(theta(loopi)) cos(theta(loopi))]; 
    t=[lmd zeros(2); zeros(2) lmd];
    dk=t'*ek*t/el(loopi);
    dm=t'*em*t*el(loopi);
    for jx=1:4
        for jy=1:4
            if(e2s(jx)*e2s(jy)~=0)
                k(e2s(jx),e2s(jy))=k(e2s(jx),e2s(jy))+dk(jx,jy);
                m(e2s(jx),e2s(jy))=m(e2s(jx),e2s(jy))+dm(jx,jy);
            end
        end
    end
end
% Inverse iterative method
p=9;
epsilon=1e-5;
[v,d]=Inviter(k,m,p,epsilon);     %d:eigenvalue  v:eigenvector  
frequency=sqrt(d)/(2*pi);
% Suspace iterative method
p=9;
epsilon=1e-5;
[v,d]=Ssiter(k,m,p,epsilon);
frequency=sqrt(d)/(2*pi);
% matlab function eigs
 [v,d] =eigs(k,m,9,'SM'); 
frequency=sqrt(diag(d))/(2*pi);
[frequency,indexf]=sort(frequency)
v=v(:,indexf);
Factor=diag(v'*m*v);
v=v*inv(sqrt(diag(Factor)));       % normalize eigenvector
