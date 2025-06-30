E=2.1e11;
A=1e-4;
density=7.3e3; 
node_number=5;
elment_number=7;
 
nc=[0,0;1,0;2,0;0,1;1,1];
%nc:node_coordinate
 
en=[1,2;2,3;1,4;2,4;2,5;3,5;4,5];
%en:element_node  
ed(1:node_number,1:2)=1;
%ed:elment_displacement 
constraint=[1,1;1,2;3,2];
 
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
 
%el:length of link element 
ek=E*A*[1 0 -1 0;0 0 0 0;-1 0 1 0;0 0 0 0];
%ek:elment_stiffness_matrix  
em=(density*A)/2*eye(4);
%em:elment_mass_matrix  
k(1:dof,1:dof)=0;
%k:structural_stiffness_matrix  
m=k;
%m:structural_mass_matrix,same size with k 
theta(1:7)=0;
el(1:7)=0;
e2s(1:4)=0;
% e2s: index of transform the elment displament number to structural
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
 
[v,d] = eig(k,m); 
%d:eigenvalue  
%v:eigenvector  
 
frequency=sqrt(diag(d))/(2*pi);

[frequency,indexf]=sort(frequency);
d=d(:,indexf);
