%-------------------------------------------------------------------------%
% Example 5.4.3
%   To solve a static deflection of a 2-D frame structure
%   nodal dofs: {u1  v1   1  u2  v2   2}
% Problem description
%   Find the deflection of a frame structure which is made of three beams
%   of lengths of 4 m, 3 m and 4 m, respectively. All beams have cross-
%   sections of 0.05 m height by 0.05 m width. The elastic modulus is
%   2.10x10^11 Pa. The frame is subjected to a concentrated load of 500 N
%   and a moment of 500 Nm at the ends of the top beam. One end of the
%   each vertical beam is fixed.
%   (see Fig. 5-9 for the element discretization)
% Variable descriptions
%   k,m = element stiffness matrix and element mass matrix
%   kk,mm = system stiffness matrix and system mass matrix
%   ff = system force vector
%   index = a vector containing system dofs associated with each element
%   bcdof = a vector containing dofs associated with boundary conditions
%   bcval = a vector containing boundary condition values associated with
%          the dofs in 'bcdof'
%   displmt = nodal displacement vector
%--------------------------------------------------------------------------
%  (0) input control data
%--------------------------------------------------------------------------
clear; clc;

Beam_InputData543;                     % import the input data for the information of
                                  % nodes, elements, loads, constraints and materials
Opt_beam=1;                                        % option for type of the beam
                                                     % =1 Euler Bernoulli beam
                                                       % =2 Timoshenko beam
Opt_mass=2;                                            % option for mass matrix
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
bcval=zeros(Sys_dof,1);                                % initializing the vector bcval

index=zeros(No_nel*No_dof,1);                         % initialization of index vector
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
for iel=1:No_el                              % loop for the total number of elements
  nd(1)=nodes(iel,1);                     % 1st connected node for the (iel)-th element
  nd(2)=nodes(iel,2);                    % 2nd connected node for the (iel)-th element
  x(1)=gcoord(nd(1),1); y(1)=gcoord(nd(1),2); z(1)=gcoord(nd(1),3);
                                                       % coordinate of 1st node
  x(2)=gcoord(nd(2),1); y(2)=gcoord(nd(2),2); z(2)=gcoord(nd(2),3);
                                                      % coordinate of 2nd node
  leng=sqrt((x(2)-x(1))^2+(y(2)-y(1))^2+(z(2)-z(1))^2);
                                                        % length of element 'iel'
                              % compute the angle between the local and global axes
  if (x(2)-x(1))==0;
    beta=pi/2; 
  elseif ((y(2)-y(1))==0)&((x(2)-x(1))<0)
    beta=pi;
  else
    beta=atan((y(2)-y(1))/(x(2)-x(1)));
  end 

  if Opt_beam==1
    [k,m]=FrameElement21(prop,leng,beta,Opt_section,Opt_mass);
                                 % compute element matrix for Euler-Bernoulli beam
  elseif Opt_beam==2
    [k,m]=FrameElement22(prop,leng,beta,Opt_section,Opt_mass);
                                   % compute element matrix for Timoshenko beam
  end

  index=femEldof(nd,No_nel,No_dof);
                                      % extract system dofs associated with element
  kk=femAssemble1(kk,k,index);                   % assemble system stiffness matrix
  mm=femAssemble1(mm,m,index);                   % assemble system mass matrix
end
%--------------------------------------------------------------------------
%  (5) apply constraints and solve the matrix equation
%--------------------------------------------------------------------------
[kk1,mm1,ff1,sdof1]=femApplybc1(kk,mm,ff,bcdof,bcval);
                             % apply fixed boundary conditions to the system equation
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
    nd(1)=nodes(iel,1);                   % 1st connected node for the (iel)-th element
    nd(2)=nodes(iel,2);                  % 2nd connected node for the (iel)-th element
    x(1)=gcoord(nd(1),1); y(1)=gcoord(nd(1),2); z(1)=gcoord(nd(1),3);
                                                       % coordinate of 1st node
    x(2)=gcoord(nd(2),1); y(2)=gcoord(nd(2),2); z(2)=gcoord(nd(2),3);
                                                      % coordinate of 2nd node
    figure(1)
      plot(x,y), xlabel('x'), ylabel('y'), hold on;
      axis([-1.5,4.5,-1,5]);
  end

  if Opt_graphics2~=1
    title('nodal connectivity of elements'), hold off;
  end
end
%---------------------------------
%  (6.2) display the static deformation
%---------------------------------
if Opt_graphics2==1
  Ampl=10;
  gcoordA=gcoord+[Ampl*displmtnode(:,1:2),zeros(No_node,1)];
  for iel=1:No_el                             % loop for the total number of elements
    nd(1)=nodes(iel,1);                   % 1st connected node for the (iel)-th element
    nd(2)=nodes(iel,2);                  % 2nd connected node for the (iel)-th element
    x(1)=gcoordA(nd(1),1); y(1)=gcoordA(nd(1),2); z(1)=gcoordA(nd(1),3);
                                                       % coordinate of 1st node
    x(2)=gcoordA(nd(2),1); y(2)=gcoordA(nd(2),2); z(2)=gcoordA(nd(2),3);
                                                      % coordinate of 2nd node
    figure(1)
      plot(x,y,'r'), title('Deflection of the structure'), hold on;
  end
  hold off
end
%--------------------------------------------------------------------------
%  (7) print fem solutions
%--------------------------------------------------------------------------
disp('The calculation is use of')

if Opt_beam==1
  disp('Euler-Bernoulli beam element')
elseif Opt_beam==2
  disp('Timoshenko beam element')
end
disp(' ')
disp('                  numrical solution')
disp('     node        x        y          ')
num=1:1:No_node;
displmts=[num' displmtnode]                             % print nodal displacements
%--------------------------------------------------------------------------
%     The end
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% Input data for static analysis of frame structure (5.4.3)
% -------------------------------------------------------------------------
%  (0.1) input basic data
% -------------------------------------------------------------------------
Prob_title='static analysis of a 2-D frame structure';

No_el=11;                                                % number of elements
No_nel=2;                                         % number of nodes per element
No_dof=3;                                            % number of dofs per node
No_node=12;                                    % total number of nodes in system
Sys_dof=No_node*No_dof;                                     % total system dofs
%--------------------------------------------------------------------------
%  (0.2) material and geometric properties
%--------------------------------------------------------------------------
prop(1)=2.1e11;                                               % elastic modulus
prop(2)=0.3;                                                   % Poisson's ratio
prop(3)=7860;                                % mass density (mass per init volume)
prop(4)=0.05;                           % height of beam cross-section in y direction
prop(5)=0.05;                            % width of beam cross-section in z direction
prop(6)=0.05*0.05;                               % cross-sectional area of the beam
prop(7)=0.05*0.05^3/12;               % moment of inertia of cross-section about axis y
prop(8)=0.05*0.05^3/12;               % moment of inertia of cross-section about axis z
 
prop(9)=1;                                             % polar moment of inertia
                   % = 0 not include the torsion deformation
                   % = 1 polar moment of inertia calculated by h*b^3 or pi*(D^4-d^4)/32
                 % = a value for the torsion deformation
prop(10)=1;                                                    % shear modulus
                                      % = 0 not include the shear deformation
                                      % = 1 shear modulus calculated by E/(2*(1+u))
                                      % = the value for the shear deformation
Opt_section=1;                          % option for correction factor for shear energy
                                                  % = 1 rectangular cross-section
                                                     % = 2 circular cross-section
%--------------------------------------------------------------------------
%  (0.3) nodal coordinates
%--------------------------------------------------------------------------
                                    % x, y, z coordinates in global coordinate system
gcoord=[ 1    0.0      0.0      0.0;
        2    0.0      1.0      0.0;
        3    0.0      2.0      0.0;
        4    0.0      3.0      0.0;
        5    0.0      4.0      0.0;
        6    1.0      4.0      0.0;
        7    2.0      4.0      0.0;
        8    3.0      4.0      0.0;
        9    3.0      3.0      0.0;
       10    3.0      2.0      0.0;
       11    3.0      1.0      0.0;
       12    3.0      0.0      0.0];
gcoord=[gcoord(:,2:4)];
%--------------------------------------------------------------------------
%  (0.4) nodal connectivity of the elements
%--------------------------------------------------------------------------
nodes=[  1          1          2;
        2          2          3;
        3          3          4;
        4          4          5;
        5          5          6;
        6          6          7;
        7          7          8;
        8          8          9;
        9          9         10;
       10         10         11;
       11         11         12];
nodes=[nodes(:,2:3)];
%--------------------------------------------------------------------------
%  (0.5) applied loads
%--------------------------------------------------------------------------
                       % a load applied at node 5 in -x direction and a moment at node 8
P=[-500,5,1;...
   -500,8,3];
%--------------------------------------------------------------------------
%  (0.6) boundary conditions
%--------------------------------------------------------------------------
ConNode=[ 1,  1,  1,  1;...                            % code for constrained nodes
         12,  1,  1,  1];
ConVal =[ 1,  0,  0,  0;...                              % values at constrained dofs
         12,  0,  0,  0];
% -------------------------------------------------------------------------
%     The end
% -------------------------------------------------------------------------
