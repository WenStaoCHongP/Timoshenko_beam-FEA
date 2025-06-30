%--------------------------------------------------------------------------
% Example 5.4.1 
%--------------------------------------------------------------------------
%   To solve the static response for a beam structure 
%   using Euler Bernoulli beam or Timoshenko beam element. 
% Problem description 
%   Find the deflection of a simply supported beam whose length is 1 m. The beam has 
%   elastic modulus of 2.1x10^11 and its cross-section is 0.02 m height by 0.02 m width. 
%   The mass density is 7860 kg/m^3. A concentrated load of -500 N is applied at the 
%   middle of the beam. 
% Variable descriptions 
%   k, m - element stiffness matrix and element mass matrix 
%   kk, mm - system stiffness matrix and system mass matrix 
%   ff - system force vector 
%   index - a vector containing system dofs associated with each element 
%   bcdof - a vector containing dofs associated with boundary conditions 
%   bcval - a vector containing boundary condition values associated with the dofs in bcdof 
%--------------------------------------------------------------------------
%  (0) input data
%--------------------------------------------------------------------------
clear; clc;

Beam_InputData541;                     % import the input data for the information of
                                  % nodes, elements, loads, constraints and materials
Opt_beam=1;                           % option for type of the beam:
                                                     % =1 Euler Bernoulli beam
                                                       % =2 Timoshenko beam
Opt_mass=2;                           % option for mass matrix:
                                                    % =1 consistent mass matrix
                                                      % =2 lumped mass matrix
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

index=zeros(No_nel*No_dof,1);                        % initialization of index vector

bcdof=zeros(Sys_dof,1);                               % initializing the vector bcdof
bcval=zeros(Sys_dof,1);                                % initializing the vector bcval
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
%  (3) applied nodal forces
%--------------------------------------------------------------------------

ff(No_dof*(P(2)-1)+P(3))=P(1);

%--------------------------------------------------------------------------
%  (4) calculate the element matrices and assembling
%--------------------------------------------------------------------------
for iel=1:No_el                               % loop for the total number of elements
  nd(1)=iel;                                 % 1st node number of the iel-th element
  nd(2)=iel+1;                              % 2nd node number of the iel-th element
  x(1)=gcoord(nd(1),1); y(1)=gcoord(nd(1),2); z(1)=gcoord(nd(1),3);
                                                       % coordinates of 1st node
  x(2)=gcoord(nd(2),1); y(2)=gcoord(nd(2),2); z(2)=gcoord(nd(2),3);
                                                      % coordinates of 2nd node
  leng=sqrt((x(2)-x(1))^2+(y(2)-y(1))^2+(z(2)-z(1))^2);
                                                        % length of element 'iel'
  if Opt_beam==1
    [k,m]=BeamElement1(prop,leng,Opt_mass);
                               % compute element matrices for Euler-Bernoulli beam
  else
    [k,m]=BeamElement2(prop,leng,Opt_mass);
                                 % compute element matrices for Timoshenko beam
  end

  index=femEldof(nd,No_nel,No_dof);      % extract system dofs associated with element

  kk=femAssemble1(kk,k,index);                   % assemble system stiffness matrix
  mm=femAssemble1(mm,m,index);                % assemble system mass matrix
end
%--------------------------------------------------------------------------
%  (5) apply constraints and solve the matrix equation
%--------------------------------------------------------------------------
[kk1,mm1,ff1]=femApplybc1(kk,mm,ff,bcdof,bcval);
                                 % apply boundary conditions to the system equation
displmt=kk1\ff1;
                                                     % solve the matrix equation
di=Sys_dof/No_dof;
for ii=1:di
  for jj=1:No_dof
    displmtnode(ii,jj)=displmt(No_dof*(ii-1)+jj,1);
  end
end
%--------------------------------------------------------------------------
%  (6) Analytical solution
%--------------------------------------------------------------------------
E=prop(1); Iz=prop(8); L=1;
nk=ceil(No_node/2);

for ii=1:nk
  Lx=(ii-1)*dx;
  c=P(1)/(48*E*Iz);
  Asolution(ii,1)=c*(3*L^2-4*Lx^2)*Lx;
  Asolution(ii,2)=c*(3*L^2-12*Lx^2);
end
for ii=nk+1:No_node
  Lx=(ii-1)*dx;
  c=P(1)/(48*E*Iz);
  Asolution(ii,1)=c*((3*L^2-4*Lx^2)*Lx+(2*Lx-L)^3);
  Asolution(ii,2)=c*((3*L^2-12*Lx^2)+6*(2*Lx-L)^2);
end
%--------------------------------------------------------------------------
%  (7) graphics of nodal connectivity and static deformation
%--------------------------------------------------------------------------
%---------------------------------
%  (7.1) display the nodal connectivity
%---------------------------------
if Opt_graphics1==1
  for iel=1:No_el                             % loop for the total number of elements
    nd(1)=iel;                               % 1st node number of the iel-th element
    nd(2)=iel+1;                            % 2nd node number of the iel-th element
    x(1)=gcoord(nd(1),1); y(1)=gcoord(nd(1),2); z(1)=gcoord(nd(1),3);
                                                       % coordinates of 1st node
    x(2)=gcoord(nd(2),1); y(2)=gcoord(nd(2),2); z(2)=gcoord(nd(2),3);
                                                      % coordinates of 2nd node
    figure(1)
      plot(x,y), xlabel('x'), ylabel('y'), hold on;
      axis([-0.1,1.1,-0.01,0.01]);
  end
  if Opt_graphics2~=1
    title('nodal connectivity of elements'), hold off;
  end
end
%---------------------------------
%  (7.2) display the static deformation
%---------------------------------
if Opt_graphics2==1
  gcoordA=gcoord+[zeros(No_node,1),displmtnode(:,1),zeros(No_node,1)];
  for iel=1:No_el                             % loop for the total number of elements
    nd(1)=iel;                               % 1st node number of the iel-th element
    nd(2)=iel+1;                            % 2nd node number of the iel-th element
    x(1)=gcoordA(nd(1),1); y(1)=gcoordA(nd(1),2); z(1)=gcoordA(nd(1),3);
                                                      % coordinates of 1st node
    x(2)=gcoordA(nd(2),1); y(2)=gcoordA(nd(2),2); z(2)=gcoordA(nd(2),3);
                                                      % coordinates of 2nd node
    figure(1)
      plot(x,y,'r'), title('Deflection of the structure'), hold on;
  end
  hold off
end
%--------------------------------------------------------------------------
%  (8) print fem solutions
%--------------------------------------------------------------------------
disp('The calculation is use of:')
if Opt_beam==1
  disp('Euler-Bernoulli beam')
else
  disp('Timoshenko beam element')
end
disp(' ')
disp('                  numrical          analytical')
disp('     node        y      thera          y       thera      ')
num=1:di;
displacements=[num' displmtnode Asolution]                % print nodal displacements
%--------------------------------------------------------------------------
%     The end
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% Input data for static analysis of beam structure (Example 5.4.1)
% -------------------------------------------------------------------------
%  (0.1) input basic data
% -------------------------------------------------------------------------
Prob_title='static analysis of beam structure';

No_el=40;                                                 % number of elements
No_nel=2;                                         % number of nodes per element
No_dof=2;                                             % number of dofs per node
No_node=(No_nel-1)*No_el+1;                      % total number of nodes in system
Sys_dof=No_node*No_dof;                                     % total system dofs
%--------------------------------------------------------------------------
%  (0.2) physical and geometric properties
%--------------------------------------------------------------------------
prop(1)=2.1e11;                                               % elastic modulus
prop(2)=0.3;                                                   % Poisson's ratio
prop(3)=7860;                                % mass density (mass per unit volume)
prop(4)=0.02;                            % height of beam cross-section in y direction
prop(5)=0.02;                             % width of beam cross-section in z direction
prop(6)=0.02*0.02;                                % cross-sectional area of the beam
prop(7)=0.02*0.02^3/12;                % moment of inertia of cross-section about axis y
prop(8)=0.02*0.02^3/12;                % moment of inertia of cross-section about axis z
%--------------------------------------------------------------------------
%  (0.3) nodal coordinates
%--------------------------------------------------------------------------
                                    % x, y, z coordinates in global coordinate system
xx=zeros(No_node,1); yy=zeros(No_node,1); zz=zeros(No_node,1);
dx=0.025; 
xx=0:dx:(No_node-1)*dx; xx=xx';
gcoord=[xx  yy  zz];
%--------------------------------------------------------------------------
%  (0.4) applied load
%--------------------------------------------------------------------------

P=[-500,21,1];                        % load applied at node 21 in the negative y direction

%--------------------------------------------------------------------------
%  (0.5) boundary conditions
%--------------------------------------------------------------------------
ConNode=[ 1, 1, 0;...                                  % code for constrained nodes
         41, 1, 0];
ConVal =[ 1, 0, 0;...                                    % values at constrained dofs
         41, 0, 0];
%--------------------------------------------------------------------------
%     The end
% -------------------------------------------------------------------------
