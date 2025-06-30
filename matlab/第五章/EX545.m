%--------------------------------------------------------------------------
% Example 5.4.5
%   To solve the eigenvalue problem for a beam structure.
%   using Euler Bernoulli beam or Timoshenko beam elements.
%   nodal dofs: {v1     v2   2}
% Problem description
%   Find the natural frequencies and modal shapes of a simply supported beam whose length
%   is 1 m. The beam has elastic modulus of 2.1x10^11 and its cross-section is 0.02 m height
%   by 0.02 m width. The mass density is 7860 kg/m^3. A concentrated load of -500 N is
%   applied at the middle of the beam.
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
Opt_beam=1;                                       % option for type of the beam:
                                                    % =1 Euler Bernoulli beam
                                                    % =2 Timoshenko beam
Opt_mass=1;                                            % option for mass matrix:
                                                    % =1 consistent mass matrix
                                                    % =2 lumped mass matrix
Opt_graphics1=0;                       % option for graphics of the nodal connectivity
Opt_graphics3=1;                            % option for graphics of the modal shape

ith=3;                                        % select the order of modes to display
%--------------------------------------------------------------------------
%  (1) initialization of matrices and vectors to zero
%--------------------------------------------------------------------------
k=zeros(No_nel*No_dof,No_nel*No_dof);                   % element stiffness matrix
m=zeros(No_nel*No_dof,No_nel*No_dof);                     % element mass matrix
 
kk=zeros(Sys_dof,Sys_dof);                   % initialization of system stiffness matrix
mm=zeros(Sys_dof,Sys_dof);                     % initialization of system mass matrix
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
%  (3) applied nodal loads
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
    [k,m]=BeamElement11(prop,leng,Opt_mass);
                               % compute element matrices for Euler-Bernoulli beam
  else
    [k,m]=BeamElement12(prop,leng,Opt_mass);
                                  % compute element matrices for Timoshenko beam
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
                     % check the boundary conditions in equation and eliminate the rows
                              % and columns associated with the boundary conditions
[V,D]=eig(kk2,mm2);                        % solve the eigenvalue problem of matrix
[lambda,ki]=sort(diag(D));                      % sort the eigenvaules and eigenvectors
omega=sqrt(lambda);                             % the frequency vector in radon/sec.
omega1=sqrt(lambda)/(2*pi);                            % the frequency vector in Hz.
V=V(:,ki);                                                   % the modal matrix
%--------------------------------------------------------------------------
%  (6) Analytical solution
%--------------------------------------------------------------------------
E=prop(1); Iz=prop(8); rho=prop(3)*prop(6); L=1;
i=(1:sdof2)';
omega2=i.*i*pi^2*sqrt(E*Iz/(rho*L^4));
omega3=omega2/(2*pi);
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
      plot(x,y); xlabel('x'), ylabel('y'), hold on;
      axis([-0.1,1.1,-1.5,1.5]);
      title('nodal connectivity of elements');
  end
  hold off
end
%--------------------------------------------------------------------------
%  (7.3) draw the modal shape
%--------------------------------------------------------------------------
if Opt_graphics3==1
  jk=ith; Vi=[V(:,jk)];                          % chose the order of modes to display
  ik=0;                        % initial the counter for locating the modal coordinates

  for ii=1:sdof0                                 % loop for the total number of sdof
    if bcdof0(ii)==0
      mcoord(ii,1)=Vi(ii-ik); 
    else
      mcoord(ii,1)=0;
      ik=ik+1;
    end
  end

  for ii=1:No_node                              % loop for the total number of nodes
    for ij=1:No_dof1
      mcoordA(ii,ij)=mcoord((ii-1)*No_dof1+ij,1);
    end
  end

  nv=20;

  for i=1:nv+1
    t=(i-1)*(2*pi)/20;
    mcoordB=gcoord+[zeros(No_node,1),mcoordA(:,1),zeros(No_node,1)]*cos(t);
    for iel=1:No_el                            % loop for the total number of elements
      nd(1)=iel;                              % starting node number for element 'iel'
      nd(2)=iel+1;                            % ending node number for element 'iel'
      x(1)=mcoordB(nd(1),1); y(1)=mcoordB(nd(1),2); z(1)=mcoordB(nd(1),3);
                                                        % coordinate of 1st node
      x(2)=mcoordB(nd(2),1); y(2)=mcoordB(nd(2),2); z(2)=mcoordB(nd(2),3);
                                                       % coordinate of 2nd node
      figure(2)
        plot(x,y,'b'), xlabel('x'), ylabel('y'), hold on;
        axis([-0.1,1.1,-1.5,1.5])
        title([num2str(jk), 'th modal shape of the structure']);
    end
    hold off
    M(:,i)=getframe;
  end
movie(M,5,10);
end
%--------------------------------------------------------------------------
%  (8) print fem solutions
%--------------------------------------------------------------------------
disp('The calculation is use of:')

if Opt_beam==1
  disp('Euler-Bernoulli beam element')
else
  disp('Timoshenko beam element')
end

if Opt_mass==1
  disp('and consistent mass matrix')
elseif Opt_mass==2
  disp('and lumped mass matrix')
end

disp(' ')
disp('     mode   numrical   analytical')

num=1:1:10;                                           % print natural frequencies
frequency=[num' omega1(1:10)  omega3(1:10)]
%--------------------------------------------------------------------------
%     The end
% -------------------------------------------------------------------------
