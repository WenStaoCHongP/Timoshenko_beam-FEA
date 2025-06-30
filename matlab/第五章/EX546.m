%--------------------------------------------------------------------------
% Example 5.4.6
%   To solve the dynamic response of a space frame structure.
%   nodal dof: {u1 v1 w1  x1  y1  z1 u2 v2 w2  x2  y2  z2}
% Problem description
%   Find the dynamic response of a space frame structure as shown in Fig.5-10, which is
%   made of twelve beam members. The frame is of the length of 4 m, the width of 3 m and
%   the height of 2 m. All members have cross-section of 0.01 m height by 0.01 m width.
%   The elastic modulus of the structure material is 2.1e11 Pa. The mass density is 7860
%   kg/m3. The frame is subjected to a concentrated load of 10000N at node 7 in -y direction.
% Variable descriptions
%   k,m = element stiffness matrix and element mass matrix
%   kk,mm = system stiffness matrix and system mass matrix
%   ff = system force vector
%   index = a vector containing system dofs associated with each element
%   bcdof = a vector containing dofs associated with boundary conditions
%   bcval = a vector containing boundary condition values associated with the dofs in 'bcdof'                                           
%--------------------------------------------------------------------------
%  (0) input control data
%--------------------------------------------------------------------------
clear; clc;

Beam_InputData544;                     % import the input data for the information of
                                  % nodes, elements, loads, constraints and materials
Opt_beam=1;                                        % option for type of the beam
                                                     % =1 Euler Bernoulli beam
                                                     % =2 Timoshenko beam
Opt_mass=1;                                            % option for mass matrix
                                                    % =1 consistent mass matrix
                                                      % =2 lumped mass matrix
Opt_section=1;                                    % option for type of cross-section
                                                  % = 1 rectangular cross-section
                                                  % = 2 circular cross-section
TypeResp=1;                           % option for selecting the type of the response
                                        % =1 impulse response:  ImpulseRespt.m
                                        % =2 step response:     StepRespt.m
                                        % =3 harmonic response: HarmonicRespt.m

dt=0.0005;                                                     % time step size
ti=0;                                                             % initial time
tf=0.25;                                                           % final time
nt=fix((tf-ti)/dt);                                           % number of time steps
tt=ti:dt:ti+nt*dt;                                            % generate time vector

ac=0.0002; bc=0.00008;                        % Parameters for proportional damping

al=0;                            % angle between the reference coordinate system and
                                  % the local coordinate system for the space element
omega0=100;                                      % frequency of the applied force
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
  xi(1)=(x(2)-x(1))/leng;                    % cos of the local x-axis and system x-axis
  xi(2)=(y(2)-y(1))/leng;                    % cos of the local x-axis and system y-axis
  xi(3)=(z(2)-z(1))/leng;                    % cos of the local x-axis and system y-axis

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
                                             % check the zero main elements in kk
[kk1,mm1,ff1,sdof1]=femApplybc1(kk0,mm0,ff0,bcdof0,bcval0);
                                 % apply boundary conditions to the system equation
 [kk2,mm2,ff2,sdof2]=bcCheck1(kk0,mm0,ff0,bcdof0,bcval0);
                         % check the boundary conditions in equation and eliminate the
                          % rows and columns associated with the boundary conditions
%--------------------------------------------------------------------------
[V,D]=eig(kk2,mm2);                        % solve the eigenvalue problem of matrix
[lambda,ki]=sort(diag(D));                     % sort the eigenvaules and eigenvectors
omega=sqrt(lambda);                                 % the frequency vector in rad/s
omega1=sqrt(lambda)/(2*pi);                            % the frequency vector in Hz

V=V(:,ki);                                                   % the modal matrix
%--------------------------------------------------------------------------
% (6) calculate dynamic response
%--------------------------------------------------------------------------
                                        % present the initial computation conditions
u=1;
C=eye(sdof2);
q0=zeros(sdof2,1);                                          % initial displacement
dq0=zeros(sdof2,1);                                             % initial velocity

if TypeResp==1
  [eta,y,omega2,sdof3]=ImpulseRespt(kk2,mm2,ff2,u,tt,C,q0,dq0,ac,bc);
                                                % compute the impulse responses
elseif TypeResp==2
  [eta,y,omega2,sdof3]=StepRespt(kk2,mm2,ff2,u,tt,C,q0,dq0,ac,bc);
                                                   % compute the step responses
else
  [eta,y,omega2,sdof3]=HarmonicRespt(kk2,mm2,ff2,omega0,tt,C,q0,dq0,ac,bc);
                                               % compute the harmonic responses
end
%--------------------------------------------------------------------------
%  (7) graphics of dynamic response
%--------------------------------------------------------------------------
jth=14;  tpoint=2;
                   % Plot the graph of the time-history response at jth degree of freedom
plot(tt,y(jth,:))
xlabel('Time (seconds)')
ylabel('Nodal  displmt. (m)')
%--------------------------------------------------------------------------
%  (8) print fem solutions
%--------------------------------------------------------------------------
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
frequency=[num' omega1]                                 % print natural frequency
%--------------------------------------------------------------------------
%     The end
% -------------------------------------------------------------------------
