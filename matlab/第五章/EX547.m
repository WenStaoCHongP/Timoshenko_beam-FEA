%-------------------------------------------------------------------------%
% Example 5.4.7
%   To solve transient response of a 2-d frame structure.
%   The solution methods are: 1) central difference scheme. 3) Houbolt integration scheme.
%   4) Wilson ?? integration scheme. 5) Newmark integration scheme
%   nodal dof: {u1  v1  w1  x1  y1  z1  u2  v2  w2  x2  y2  z2}
% Problem description
%   Find the response of a frame structure which is made of three beams of lengths of 4 m, 
%   3 m and 4 m, respectively. All beams have cross- section of 0.10 m height by 0.05 m
%   width. The elastic modulus is 2.10x10^11 Pa. The frame is subjected to an impulse load
%   of amplitude 500 N in the middle of the top beam. One end of the each vertical beam is
%   fixed. (see Fig. 5-9 for the element discretization)
% Variable descriptions
%   k, m - element stiffness matrix and mass matrix
%   kk, mm - system stiffness matrix and mass matrix
%   ff - system force vector
%   index - a vector containing system dofs associated with each element
%   bcdof - a vector containing dofs associated with boundary conditions
%   bcval - a vector containing boundary condition values associated with the dofs in 'bcdof'
%   dsp - displacement matrix
%   vel - velocity matrix
%   acc - acceleartion matrix
%--------------------------------------------------------------------------
%  (0) input control data
%--------------------------------------------------------------------------
clear; clc;

Beam_InputData547;                     % import the input data for the information of
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
TypeMethod=1;                            % option for selecting the solution method
                                                 % = 1 central difference scheme
                                               % = 3 Houbolt integration scheme
                                             % = 4 Wilson   integration scheme
                                              % = 5 Newmark integration scheme
Typeload=1;                                    % option for selecting the load type
                                                           % = 1 impulse load
                                                              % = 2 step load
                                                          % = 3 Harmonic load

dt=0.00001;                                                    % time step size
ti=0;                                                            % initial time
tf=0.200;                                                          % final time
nt=fix((tf-ti)/dt);                                           % number of time steps
tt=ti:dt:ti+nt*dt;                                     % generate time samples vector
 
ac=0.00002; bc=0.00008;                       % Parameters for proportional damping
 
al=0;                            % angle between the reference coordinate system and
                                  % the local coordinate system for the space element
%--------------------------------------------------------------------------
%  (1) initialization of matrices and vectors to zero
%--------------------------------------------------------------------------
k=zeros(No_nel*No_dof,No_nel*No_dof);                   % element stiffness matrix
m=zeros(No_nel*No_dof,No_nel*No_dof);                     % element mass matrix
 
kk=zeros(Sys_dof,Sys_dof);                   % initialization of system stiffness matrix
mm=zeros(Sys_dof,Sys_dof);                     % initialization of system mass matrix
ff=zeros(Sys_dof,1);                            % initialization of system force vector
 
bcdof=zeros(Sys_dof,1);                               % initializing the vector bcdof
bcval=zeros(Sys_dof,1);                                % initializing the vector bcval
 
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
  nd(1)=iel;                                 % starting node number for element 'iel'
  nd(2)=iel+1;                                % ending node number for element 'iel'
  x(1)=gcoord(nd(1),1); y(1)=gcoord(nd(1),2); z(1)=gcoord(nd(1),3);
                                                       % coordinate of 1st node
  x(2)=gcoord(nd(2),1); y(2)=gcoord(nd(2),2); z(2)=gcoord(nd(2),3);
                                                      % coordinate of 2nd node
  leng=sqrt((x(2)-x(1))^2+(y(2)-y(1))^2+(z(2)-z(1))^2);
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

  index=femEldof(nd,No_nel,No_dof);      % extract system dofs associated with element
  kk=femAssemble1(kk,k,index);                   % assemble system stiffness matrix
  mm=femAssemble1(mm,m,index);                   % assemble system mass matrix
 
end
%--------------------------------------------------------------------------
%  (5) apply constraints and solve the matrix equation
%--------------------------------------------------------------------------
[kk0,mm0,ff0,bcdof0,bcval0,sdof0]=kkCheck1(kk,mm,ff,bcdof,bcval);
                          % check the zero main elements in kk and eliminate the rows
                      % and columns in equation associated with the zero main elements
[kk1,mm1,ff1,sdof1]=femApplybc1(kk0,mm0,ff0,bcdof0,bcval0);
                                 % apply boundary conditions to the system equation

[kk2,mm2,ff2,sdof2]=bcCheck1(kk0,mm0,ff0,bcdof0,bcval0);
                            % check the boundary conditions in equation and eliminate
                       % the rows and columns associated with the boundary conditions

[V,D]=eig(kk2,mm2);                        % solve the eigenvalue problem of matrix
[lambda,ki]=sort(diag(D));                       % sort the eigenvaules and eigenvectors
omega=sqrt(lambda);                                % the frequency vector in radin/s
omega1=sqrt(lambda)/(2*pi);                             % the frequency vector in Hz

V=V(:,ki);                                                   % the modal matrix
%--------------------------------------------------------------------------
%  (6) Check the eigenvalues and the damping parameters
%--------------------------------------------------------------------------
                            % check whether the eigenvalues are infinite and eliminate
                               % the eigenvectors associated with the bad eigenvalues.
jk=0;

for ii=1:sdof2                                   % loop for find the infinite in omega
  check=omega(ii);
  if check>1.0e12
    jk=jk+1;                                    % location of the infinite frequency
    omi(jk)=ii;                         % storing the location of the infinite frequency
  end
end

sdof3=sdof2-jk;
V1=[V(:,1:sdof3)];                                     % truncate the modal vectors

Factor=diag(V1'*mm2*V1);
Vnorm=V1*inv(sqrt(diag(Factor)));                    %  Eigenvectors are normalized
%--------------------------------------------------------------------------
omega2=diag(sqrt(Vnorm'*kk2*Vnorm));                       % Natural frequencies

Modamp=Vnorm'*(ac*mm2+bc*kk2)*Vnorm;             % Form the Rayleigh damping
zeta=diag((1/2)*Modamp*inv(diag(omega2)));                     % The damping ratio

if (max(zeta) >= 1)
  disp('Warning - Your maximu damping ratio is grater than or equal to 1')
  disp('You have to reselect ac and bc ')
  pause
  disp('If you want to continue, type return key')
end
%--------------------------------------------------------------------------
%  (7) calculate transient response
%--------------------------------------------------------------------------
[kk1,mm1,ff1,bcdof1,bcval1,sdof1]=mmCheck1(kk0,mm0,ff0,bcdof0,bcval0);
                         % check the zero main elements in mm and eliminate the rows
                      % and columns in equation associated with the zero main elements
switch Typeload
  case 1                                                % Impulse force function
    u=[1,zeros(1,nt)];
    ft0=ff1*u;
  case 2                                                   % Step force function
    u(1,1:nt+1)=1;
    ft0=ff1*u;
  case 3                                              % Harmonic force function
    u=cos(omega0*tt);
    ft0=ff1*u;
  otherwise
    ft0=ff1;                                             % a given force function
end

cc1=ac*mm1+bc*kk1;                          % Form the Rayleigh damping matrix

q0=zeros(sdof1,1); dq0=zeros(sdof1,1);  % initial displacement and velocity
 
switch TypeMethod
  case 1                                            % central difference scheme 1
    [acc,vel,dsp]=TransResp1(kk1,cc1,mm1,ft0,bcdof1,nt,dt,q0,dq0);
  case 3                                            % Houbolt integration scheme
    [acc,vel,dsp]=TransResp3(kk1,cc1,mm1,ft0,bcdof1,nt,dt,q0,dq0);
  case 4                                           % Wilson ?? integration scheme
    [acc,vel,dsp]=TransResp4(kk1,cc1,mm1,ft0,bcdof1,nt,dt,q0,dq0);
  case 5                                           % Newmark integration scheme
    [acc,vel,dsp]=TransResp5(kk1,cc1,mm1,ft0,bcdof1,nt,dt,q0,dq0);
  otherwise
    disp('Unknown method.')
end
%--------------------------------------------------------------------------
%  (8) graphics of dynamic response
%--------------------------------------------------------------------------
  jth=14;
                   % Plot the graph of the time-history response at jth degree of freedom
  plot(tt,dsp(jth,:))
  xlabel('Time  (seconds)'), ylabel('displacement  (m)')
  title('time-history response')
%--------------------------------------------------------------------------
%  (9) print fem solutions
%--------------------------------------------------------------------------
switch Typeload
  case 1
    disp('The excitation is impulse force')
  case 2
    disp('The excitation is step force')
  case 3
    disp('The excitation is step force')
    otherwise
    disp('The given foece')
end

disp('The calculation is use of:')
 
if Opt_beam==1
  disp('Euler-Bernoulli beam element')
elseif Opt_beam==2
  disp('Timoshenko beam element')
end
 
if Opt_mass==1
  disp('and consistent mass matrix')
elseif Opt_mass==2
  disp('and lumped mass matrix')
end
 
num=1:1:sdof2;
frequency=[num' omega1]                                  % print natural frequency
%--------------------------------------------------------------------------
%     The end
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Input data for static analysis of frame structure (5.4.7)
%--------------------------------------------------------------------------
%  (0.1) input basic data
%--------------------------------------------------------------------------
Prob_title='transient analysis of a 2-D frame structure';
 
No_el=11;                                                 % number of elements
No_nel=2;                                         % number of nodes per element
No_dof=6;                                             % number of dofs per node
No_node=12;                                     % total number of nodes in system
Sys_dof=No_node*No_dof;                                     % total system dofs
%--------------------------------------------------------------------------
%  (0.2) material and geometric properties
%--------------------------------------------------------------------------
prop(1)=2.1e11;                                               % elastic modulus
prop(2)=0.3;                                                   % Poisson's ratio
prop(3)=7860;                                 % mass density (mass per init volume)
prop(4)=0.10;                            % height of beam cross-section in y direction
prop(5)=0.05;                             % width of beam cross-section in z direction
prop(6)=0.05*0.10;                                % cross-sectional area of the beam
prop(7)=0;                           % moment of inertia of cross-section about axis y
prop(8)=0.05*0.10^3/12;                % moment of inertia of cross-section about axis z

prop(9)=0;                                              % polar moment of inertia
                                        % = 0 not include the torshional deformation
                                        % = 1 polar moment of inertia calculated by
                                                %     h*b^3 or pi*(D^4-d^4)/32
                                        % = the value for the torshional deformation
prop(10)=1;                                                    % shear modulus
                                           % = 0 not include the shear deformation
                                     % = 1 shear modulus calculated by E/(2*(1+u))
                                             % = a value for the shear deformation
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
                                          % a load applied at node 6 in -y direction
P=[-500,6,2];

%--------------------------------------------------------------------------
%  (0.6) boundary conditions
%--------------------------------------------------------------------------
ConNode=[ 1,  1,  1,  1,  1,  1,  1;                   % code for constrained nodes
         12,  1,  1,  1,  1,  1,  1];
ConVal =[ 1,  0,  0,  0,  0,  0,  0;                     % values at constrained dofs
         12,  0,  0,  0,  0,  0,  0];
%--------------------------------------------------------------------------
%     The end
%--------------------------------------------------------------------------
