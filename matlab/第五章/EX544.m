%-------------------------------------------------------------------------%
% Example 5.4.4
%   To solve a static deflection of a space frame structure.
%   nodal dof: {u1  v1  w1  _x1  _y1  _z1  u2  v2  w2  _x2  _y2  _z2}
% Problem description
%   Find the static deflection of a space frame structure as shown in Fig.5-11, which is made
%   of twelve beam members. The frame are of the length of 4 m, the width of 3 m and the
%   height of 2 m. All members have cross-section of 0.01 m height by 0.01 m width. The
%   elastic modulus of the structure material is 2.1e11 Pa. The mass density is 7860 kg/m3.
%   The frame is subjected to a concentrated load of 10000N at node 7 in -y direction.
% Variable descriptions
%   k,m = element stiffness matrix and element mass matrix
%   kk,mm = system stiffness matrix and system mass matrix
%   ff = system force vector
%   index = a vector containing system dofs associated with each element
%   bcdof = a vector containing dofs associated with boundary conditions
%   bcval = a vector containing boundary condition values associated with the dofs in 'bcdof'
%   displmt = nodal displacement vector
%-------------------------------------------------------------------------%
%  (0) input control data
%--------------------------------------------------------------------------
clear; clc;

Beam_InputData544;                     % import the input data for the information of
                                  % nodes, elements, loads, constraints and materials
Opt_beam=1;                                        % option for type of the beam:
                                                     % =1 Euler Bernoulli beam
                                                       % =2 Timoshenko beam
Opt_mass=2;                                            % option for mass matrix:
                                                % =1 consistent mass matrix
                                                      % =2 lumped mass matrix
Opt_section=1;                                    % option for type of cross-section
                                                  % = 1 rectangular cross-section
                                                     % = 2 circular cross-section
Opt_graphics1=1;                       % option for graphics of the nodal connectivity
Opt_graphics2=1;                       % option for graphics of the static deformation
%--------------------------------------------------------------------------
%  (1) initialization of matrices and vectors to zero
%--------------------------------------------------------------------------
k=zeros(No_nel*No_dof,No_nel*No_dof);                   % element stiffness matrix
m=zeros(No_nel*No_dof,No_nel*No_dof);                     % element mass matrix

kk=zeros(Sys_dof,Sys_dof);                   % initialization of system stiffness matrix
mm=zeros(Sys_dof,Sys_dof);                    % initialization of system mass matrix
ff=zeros(Sys_dof,1);                            % initialization of system force vector

bcdof=zeros(Sys_dof,1);                               % initializing the vector bcdof
bcval=zeros(Sys_dof,1);                               % initializing the vector bcval

index=zeros(No_nel*No_dof,1);                        % initialization of index vector
%--------------------------------------------------------------------------
%  (2) calculation of constraints
%--------------------------------------------------------------------------
[n1,n2]=size(ConNode);
                                                 % calculate the constrained dofs
for ni=1:n1
  ki=ConNode(ni,1);
  bcdof((ki-1)*No_dof+1:ki*No_dof)=ConNode(ni,2:No_dof+1);
                                                  % the code for constrained dofs
  bcval((ki-1)*No_dof+1:ki*No_dof)=ConVal(ni,2:No_dof+1);
                                                  % the value at constrained dofs
end
%--------------------------------------------------------------------------
%  (3) applied nodal loads
%--------------------------------------------------------------------------
[n1,n2]=size(P);

for ni=1:n1
  ff(No_dof*(P(ni,2)-1)+P(ni,3))=P(ni,1);
end
%--------------------------------------------------------------------------
%  (4) calculate the element matrices and assembling
%--------------------------------------------------------------------------
for iel=1:No_el                               % loop for the total number of elements
  nd(1)=nodes(iel,1);                       % 1st connected node for the iel-th element
  nd(2)=nodes(iel,2);                      % 2nd connected node for the iel-th element
  x(1)=gcoord(nd(1),1); y(1)=gcoord(nd(1),2); z(1)=gcoord(nd(1),3);
                                                        % coordinate of 1st node
  x(2)=gcoord(nd(2),1); y(2)=gcoord(nd(2),2); z(2)=gcoord(nd(2),3);
                                                       % coordinate of 2nd node
  leng=sqrt((x(2)-x(1))^2+(y(2)-y(1))^2+(z(2)-z(1))^2);
                                                         % length of element 'iel'
  xi(1)=(x(2)-x(1))/leng;                     % cos of the local x-axis and system x-axis
  xi(2)=(y(2)-y(1))/leng;                     % cos of the local x-axis and system y-axis
  xi(3)=(z(2)-z(1))/leng;                     % cos of the local x-axis and system y-axis

  if Opt_beam==1
    [k,m]=FrameElement31(prop,leng,xi,al,Opt_section,Opt_mass);
                                 % compute element matrix for Euler-Bernoulli beam
  elseif Opt_beam==2
    [k,m]=FrameElement32(prop,leng,xi,al,Opt_section,Opt_mass);
                                    % compute element matrix for Timoshenko beam
  end

  index=femEldof(nd,No_nel,No_dof);
                              % extract system dofs associated with element
  kk=femAssemble1(kk,k,index);                   % assemble system stiffness matrix
  mm=femAssemble1(mm,m,index);                % assemble system mass matrix
end
%--------------------------------------------------------------------------
%  (5) apply constraints and solve the matrix equation
%--------------------------------------------------------------------------
[kk0,mm0,ff0,bcdof0,bcval0,sdof0]=kkCheck1(kk,mm,ff,bcdof,bcval);
                          % check the zero main elements in kk and eliminate the rows
                      % and columns in equation associated with the zero main elements
[kk1,mm1,ff1,sdof1]=femApplybc1(kk0,mm0,ff0,bcdof0,bcval0);
                                  % apply boundary conditions to the system equation

displmt=kk1\ff1;                 % solve the matrix equation to find nodal displacements

No_dof1=No_dof/(Sys_dof/sdof1);
for ii=1:No_node
  for ij=1:No_dof1
    displmtnode(ii,ij)=displmt((ii-1)*No_dof1+ij,1);
  end
end
%--------------------------------------------------------------------------
%  (6) graphics of nodal connectivity and static deformation
%--------------------------------------------------------------------------
%---------------------------------
%  (6.1) display the nodal connectivity
%---------------------------------
if Opt_graphics1==1
  for iel=1:No_el                            % loop for the total number of elements
    nd(1)=nodes(iel,1);                    % 1st connected node for the iel-th element
    nd(2)=nodes(iel,2);                   % 2nd connected node for the iel-th element
    x(1)=gcoord(nd(1),1); y(1)=gcoord(nd(1),2); z(1)=gcoord(nd(1),3);
                                                      % coordinate of 1st node
    x(2)=gcoord(nd(2),1); y(2)=gcoord(nd(2),2); z(2)=gcoord(nd(2),3);
                                                     % coordinate of 2nd node
    figure(1)
      plot3(x,y,z), xlabel('x'), ylabel('y'), zlabel('z'), hold on;
  end
  if Opt_graphics2~=1
    title('nodal connectivity of elements'); hold off;
  end
end
%---------------------------------
%  (6.2) display the static deformation
%---------------------------------
if Opt_graphics2==1
  gcoordA=gcoord+10*displmtnode(:,1:3);
  for iel=1:No_el                             % loop for the total number of elements
    nd(1)=nodes(iel,1);                     % 1st connected node for the iel-th element
    nd(2)=nodes(iel,2);                    % 2nd connected node for the iel-th element
    x(1)=gcoordA(nd(1),1); y(1)=gcoordA(nd(1),2); z(1)=gcoordA(nd(1),3);
                                                       % coordinate of 1st node
    x(2)=gcoordA(nd(2),1); y(2)=gcoordA(nd(2),2); z(2)=gcoordA(nd(2),3);
                                                      % coordinate of 2nd node
    figure(1)
      plot3(x,y,z,'r'), title('Deflection of the structure'),  hold on;
  end
  hold off
end
%--------------------------------------------------------------------------
%  (7) print fem solutions
%--------------------------------------------------------------------------
disp('The calculation is use of:')

if Opt_beam==1
  disp('Euler-Bernoulli beam element')
elseif Opt_beam==2
  disp('Timoshenko beam element')
end
disp(' ')
disp('                  numrical solution')
disp('     node        x         y        z          x        y        z ')
num=1:1:No_node;
displmts=[num' displmtnode]                             % print nodal displacements
%--------------------------------------------------------------------------
%     The end
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% Input data for static analysis of frame structure (5.4.4)
% -------------------------------------------------------------------------
%  (0.1) input basic data
% -------------------------------------------------------------------------
Prob_title='static analysis of a 3-D frame structure';

No_el=12;                                                % number of elements
No_nel=2;                                         % number of nodes per element
No_dof=6;                                            % number of dofs per node
No_node=8;                                     % total number of nodes in system
Sys_dof=No_node*No_dof;                                     % total system dofs
%--------------------------------------------------------------------------
%  (0.2) material and geometric properties
%--------------------------------------------------------------------------
prop(1)=2.1e11;                                               % elastic modulus
prop(2)=0.3;                                                   % Poisson's ratio
prop(3)=7860;                                % mass density (mass per init volume)
prop(4)=0.01;                           % height of beam cross-section in y direction
prop(5)=0.01;                            % width of beam cross-section in z direction
prop(6)=0.01*0.01;                               % cross-sectional area of the beam
prop(7)=0.01*0.01^3/12;               % moment of inertia of cross-section about axis y
prop(8)=0.01*0.01^3/12;               % moment of inertia of cross-section about axis z
 
prop(9)=1;                                             % polar moment of inertia
                                       % = 0 not include the torsion deformation
                                       % = 1 polar moment of inertia calculated by
                                       %    h*b^3 or pi*(D^4-d^4)/32
                                       % = a value for the torsion deformation
prop(10)=1;                                                   % shear modulus
                                     % = 0 not include the shear deformation
                                     % = 1 shear modulus calculated by E/(2*(1+u))
                                     % = the value for the shear deformation
Opt_section=1;                          % option for correction factor for shear energy
                                                   % = 1 rectangular cross-section
                                                   % = 2 circular cross-section
%--------------------------------------------------------------------------
%  (0.3) nodal coordinates
%--------------------------------------------------------------------------
                            % x, y, z coordinates of nodes in global coordinate system
gcoord=[ 1    0,  0,  0;
        2    4,  0,  0;
        3    4,  3,  0;
0,  3,  0;
        5    0,  0,  2;
        6    4,  0,  2;
        7    4,  3,  2;
        8    0,  3,  2];
gcoord=[gcoord(:,2:4)];
al=0;                            % angle between the reference coordinate system and
                                  % the local coordinate system for the space element
%--------------------------------------------------------------------------
%  (0.4) nodal connectivity of the elements
%--------------------------------------------------------------------------
nodes=[ 1    1,  5;
       2    2,  6;
       3    3,  7;
       4    4,  8;
       5    5,  6;
6    6,  7;
       7    7,  8;
       8    8,  5;
       9    1,  6;
      10    2,  7;
      11    3,  8;
      12    4,  5];
nodes=[nodes(:,2:3)];
%--------------------------------------------------------------------------
%  (0.5) applied loads
%--------------------------------------------------------------------------
                                          % a load applied at node 7 in -y direction
P=[-100000,7,2];

%--------------------------------------------------------------------------
%  (0.6) boundary conditions
%--------------------------------------------------------------------------
ConNode=[ 1,  1,  1,  1,  1,  1,  1;                   % code for constrained nodes
          2,  1,  1,  1,  1,  1,  1;
          3,  1,  1,  1,  1,  1,  1;
          4,  1,  1,  1,  1,  1,  1];
ConVal =[  1,  0,  0,  0,  0,  0,  0;                   % values at constrained dofs
          2,  0,  0,  0,  0,  0,  0;
          3,  0,  0,  0,  0,  0,  0;
          4,  0,  0,  0,  0,  0,  0];
% -------------------------------------------------------------------------
%     The end
% -------------------------------------------------------------------------
