function [W,dispmax,stressmax,frequency]=Trusssizeopt(area)
%--------------------------------------------------------------------------
% Purpose:
% This function is used to calculate total mass, maximal nodal diaplcement,
%  maximal element stress and first natural frequency for space 72-bar truss 
%
% Synopsis:
%       [W,dispmax,stressmax,frequency]=Trusssizeopt(area)
% Variable Description:
%       Input parameters
%       area - element cross-section areas 
%       Output parameters
%       W - total truss mass
%       dispmax - maximal nodal displacement 
%       stressmax - maximal element stress
%       frequency - first natural frequency
%--------------------------------------------------------------------------
E=2.1e11;                                % elastic modulus
density=7860;                            % mass density
node_number=20;
element_number=72;
No_dof=3;
No_nel=2;
Sdof=node_number*No_dof;                 % total number of dofs for ground structure
f=zeros(Sdof,1);                         % force vector
f(58)=-50000;
nc=[0 0 0;1 0 0;1 1 0;0 1 0; ...        % nodal coordinate for ground structure
    0 0 1;1 0 1;1 1 1;0 1 1; ...
    0 0 2;1 0 2;1 1 2;0 1 2; ...
     0 0 3;1 0 3;1 1 3;0 1 3; ...
     0 0 4;1 0 4;1 1 4;0 1 4];
en=zeros(72,2);                             
en(1,1)=1;en(1,2)=5;                                       %en:element_node    
en(2,1)=1;en(2,2)=6;
en(3,1)=2;en(3,2)=5;
en(4,1)=2;en(4,2)=6;
en(5,1)=2;en(5,2)=7;
en(6,1)=3;en(6,2)=6;
en(7,1)=3;en(7,2)=7;
en(8,1)=3;en(8,2)=8;
en(9,1)=4;en(9,2)=7;
en(10,1)=4;en(10,2)=8;
en(11,1)=4;en(11,2)=5;
en(12,1)=1;en(12,2)=8;
en(13,1)=5;en(13,2)=6;
en(14,1)=6;en(14,2)=7;
en(15,1)=7;en(15,2)=8;
en(16,1)=5;en(16,2)=8;
en(17,1)=6;en(17,2)=8;
en(18,1)=5;en(18,2)=7;
for loopi=1:3
    for loopj=1:18
       en(loopi*18+loopj,1)=en(loopj,1)+4*loopi;
       en(loopi*18+loopj,2)=en(loopj,2)+4*loopi;
    end
end

ed(1:node_number,1:3)=1;                            %ed:elment_displacement 
constraint=[1,1;1,2;1,3;2,1;2,2;2,3;3,1;3,2;3,3;4,1;4,2;4,3];
for loopi=1:length(constraint);
    ed(constraint(loopi,1),constraint(loopi,2))=0;
end
dof=0;
for loopi=1:node_number
     for loopj=1:3
        if ed(loopi,loopj)~=0
            dof=dof+1;
            ed(loopi,loopj)=dof;
        end
    end
end
%-----------------------------------------
% Initilization to zero
%-----------------------------------------
k(1:dof,1:dof)=0;                    %structural_stiffness_matrix  
m(1:dof,1:dof)=0;;                   %structural_mass_matrix
disp(1:dof)=0;                       % system displacement vector
eldisp(1:No_nel*No_dof)=0;           % element nodal displacement vector
elforce(1:No_nel*No_dof)=0;          % element force vector
stress(1:element_number)=0;           % stress vector for every element
el(1:element_number)=0;
e2s(1:6)=0;      % index of transform the elment displament number to structural
%---------------------------------------------------------
% Assemble system mass and siffness matrix
%---------------------------------------------------------
W=0;
for loopi=1:element_number
    for zi=1:2
        e2s((zi-1)*2+1)=ed(en(loopi,zi),1);
        e2s((zi-1)*2+2)=ed(en(loopi,zi),2);
        e2s((zi-1)*2+3)=ed(en(loopi,zi),3);
    end
    el(loopi)=sqrt((nc(en(loopi,1),1)-nc(en(loopi,2),1))^2+(nc(en(loopi,1),2)-nc(en(loopi,2),2))^2 ...
              +(nc(en(loopi,1),3)-nc(en(loopi,2),3))^2);
    c=(nc(en(loopi,1),1)-nc(en(loopi,2),1))/el(loopi);
    s=(nc(en(loopi,1),2)-nc(en(loopi,2),2))/el(loopi);
    r=(nc(en(loopi,1),3)-nc(en(loopi,2),3))/el(loopi);
    dk=E*area(loopi)/el(loopi)*[c^2 c*s c*r -c^2 -c*s -c*r; ...   % element stiffness matrix
                                c*s s^2 s*r -c*s -s^2 -s*r; ...
                                c*r s*r r^2 -c*r -s*r -r^2; ...
                                -c^2 -c*s -c*r c^2 c*s c*r; ...
                                -c*s -s^2 -s*r c*s s^2 s*r; ...
                                -c*r -s*r -r^2 c*r s*r r^2];
    dm=(density*area(loopi)*el(loopi))/2*eye(6);                 % element mass matrix
    for jx=1:6
        for jy=1:6
             if(e2s(jx)*e2s(jy)~=0)
                k(e2s(jx),e2s(jy))=k(e2s(jx),e2s(jy))+dk(jx,jy);
                m(e2s(jx),e2s(jy))=m(e2s(jx),e2s(jy))+dm(jx,jy);
             end
        end
    end
    W=W+density*area(loopi)*el(loopi);
end
%---------------------------------------
% Compute the first natural frequency
%---------------------------------------
p=1;                             % the number of natural frequencies
epsilon=1e-5;
[v,d]=Inviter(k,m,p,epsilon);    % Inverse iterative method
frequency=sqrt(d)/(2*pi); 
%---------------------------------------
% Compute system displacement
%--------------------------------------
 Number_con=-1;                % apply boundary conditions
 for loopi=1:length(constraint)   
     Number_con= Number_con+1;
     f((constraint(loopi,1)-1)*No_dof+constraint(loopi,2)-Number_con)=[];
 end
  disp=k\f;                   % solve matrix equation for nodal displacement
%-----------------------------------------
% post computation for stress calculation
%----------------------------------------
  for loopi=1:element_number
    for zi=1:2
        e2s((zi-1)*2+1)=ed(en(loopi,zi),1);
        e2s((zi-1)*2+2)=ed(en(loopi,zi),2);
        e2s((zi-1)*2+3)=ed(en(loopi,zi),3);
    end
    el(loopi)=sqrt((nc(en(loopi,1),1)-nc(en(loopi,2),1))^2+(nc(en(loopi,1),2)-nc(en(loopi,2),2))^2 ...
              +(nc(en(loopi,1),3)-nc(en(loopi,2),3))^2);
    c=(nc(en(loopi,1),1)-nc(en(loopi,2),1))/el(loopi);
    s=(nc(en(loopi,1),2)-nc(en(loopi,2),2))/el(loopi);
    r=(nc(en(loopi,1),3)-nc(en(loopi,2),3))/el(loopi);
    dk=E*area(loopi)/el(loopi)*[c^2 c*s c*r -c^2 -c*s -c*r; ...
                                c*s s^2 s*r -c*s -s^2 -s*r; ...
                                c*r s*r r^2 -c*r -s*r -r^2; ...
                                -c^2 -c*s -c*r c^2 c*s c*r; ...
                                -c*s -s^2 -s*r c*s s^2 s*r; ...
                                -c*r -s*r -r^2 c*r s*r r^2];
    for loopj=1:(No_nel*No_dof)
       if e2s(loopj)==0
          eldisp(loopj)=0;
       else
          eldisp(loopj)=disp(e2s(loopj));
       end
    end
     elforce=dk*eldisp';                            % element force vector
    if area(loopi)==0
        stress(loopi)=0;
    else
       stress(loopi)=sqrt(elforce(1)^2+elforce(2)^2+elforce(3)^2)/area(loopi); % stress
    end
    if ((nc(en(loopi,1),1)-nc(en(loopi,2),1))*elforce(4))<0      % check if tension or compression
       stress(loopi)=-stress(loopi);
    end
  end
%---------------------------------------
% Construct the taotal displacement
%---------------------------------------
 flag=1;
 iter=-1;
 displ=disp;
 while flag
     for loopi=1:length(constraint)
         iter=iter+1;
         for loopj=1:(dof+iter)
             if loopj<((constraint(loopi,1)-1)*No_dof+constraint(loopi,2))
                 dispt(loopj)=displ(loopj);
             else
                 dispt(loopj+1)=displ(loopj);
             end
         end  
         dispt((constraint(loopi,1)-1)*No_dof+constraint(loopi,2))=0;  
         displ=dispt;
     end 
     if length(displ)>=node_number*No_dof
             flag=0;
     end 
 end
%---------------------------------------------------
%  Output solutions
%---------------------------------------------------
  dispmax=max(abs(dispt'));
  stressmax=max(abs(stress'));
