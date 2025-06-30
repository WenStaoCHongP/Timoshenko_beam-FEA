%--------------------------------------------------------------------------
% Example 4.7.2
%--------------------------------------------------------------------------
%   plane stress analysis of a cantilever beam using isoparametric 
%   four-node elements. (see Fig. 4-10 for the finite element mesh)
%--------------------------------------------------------------------------
%----------------------------------
%  (0) input data for control parameters
%----------------------------------
clear; clc;

No_element=8;                 % number of elements
No_nodeEl=4;                  % number of nodes per element
No_dof=2;                     % number of dofs per node
No_nodeSys=18;                % total number of nodes in system
Sys_dof=No_nodeSys*No_dof;    % total system dofs  
El_dof=No_nodeEl*No_dof;       % degrees of freedom per element
emodule=1e6;                   % elastic modulus
poisson=0.3;                    % Poisson's ratio
nglx=2; ngly=2;                 % 2x2 Gauss quadrature
nglxy=nglx*ngly;                % number of sampling points per element
%-------------------------------
%  (1) input nodal coordinate values
%-------------------------------
       % x, y coordinates of nodes in global coordinate system
gcoord=[  0.0  0.0;    0.0  1.0;
         0.5  0.0;    0.5  1.0;
         1.0  0.0;    1.0  1.0;
         1.5  0.0;    1.5  1.0;
         2.0  0.0;    2.0  1.0;
         2.5  0.0;    2.5  1.0;
         3.0  0.0;    3.0  1.0;
         3.5  0.0;    3.5  1.0;
         4.0  0.0;    4.0  1.0];
%------------------------------------------------
%  (2) input data for nodal connectivity for each element
%------------------------------------------------
%      element     node code
%        No.   i    j    k    m
nodes=[   1    1    3    4    2;
         2     3    5    6    4;
         3     5    7    8    6;
         4     7    9   10    8;
         5     9   11   12   10;
         6    11   13   14   12;
         7    13   15   16   14;
         8    15   17   18   16];
nodes=[nodes(:,2:5)];
%-----------------------------------
%  (3) input data for boundary conditions
%-----------------------------------
bcdof=[1  2  3  4];        % first four dofs are constrained
bcval=[0  0  0  0];        % whose described values are 0 
%------------------------------------
%  ги4) initialization of matrices and vectors
%------------------------------------
ff=zeros(Sys_dof,1);         % system force vector
kk=zeros(Sys_dof,Sys_dof);   % system matrix
displmt=zeros(Sys_dof,1);    % system displacement vector
eldisp=zeros(El_dof,1);       % element displacement vector
stress=zeros(nglxy,3);        % matrix containing stress components
strain=zeros(nglxy,3);        % matrix containing strain components
index=zeros(El_dof,1);       % index vector
kinmtx2=zeros(3,El_dof);     % kinematic matrix
matmtrx=zeros(3,3);          % constitutive matrix
%-----------------
%  (5) force vector
%-----------------
ff(34)=500;                  % force applied at node 17 in y-axis
ff(36)=500;                  % force applied at node 18 in y-axis
%----------------------------------------------------------
%  (6) computation of element matrices and vectors and their assembly
%----------------------------------------------------------
[point2,weight2]=GaussPoint2(nglx,ngly);     % sampling points & weights
matmtrx=Materialiso(1,emodule,poisson);      % compute constitutive matrix

for iel=1:No_element         % loop for the total number of elements

for i=1:No_nodeEl
nd(i)=nodes(iel,i);            % extract connected node for (iel)-th element
xcoord(i)=gcoord(nd(i),1);     % extract x value of the node
ycoord(i)=gcoord(nd(i),2);     % extract y value of the node
end

k=zeros(El_dof,El_dof);       % initialization of element matrix to zero
%-------------------------------------------------
%  (7) numerical integration and  compute element matrix
%-------------------------------------------------
for intx=1:nglx
x=point2(intx,1);                     % sampling point in x-axis
wtx=weight2(intx,1);                  % weight in x-axis
for inty=1:ngly
y=point2(inty,2);                     % sampling point in y-axis
wty=weight2(inty,2) ;                 % weight in y-axis

[shape,dhdr,dhds]=FemIsoq4(x,y);       % compute shape functions and
                                   % derivatives at sampling point
jacob2=FemJacobi2(No_nodeEl,dhdr,dhds,xcoord,ycoord);  % compute Jacobian

detjacob=det(jacob2);                 % determinant of Jacobian
invjacob=inv(jacob2);                 % inverse of Jacobian matrix

[dhdx,dhdy]=FemDeriv2(No_nodeEl,dhdr,dhds,invjacob);  % derivatives w.r.t.
                                                 % physical coordinate
kinmtx2=FemKine2D(No_nodeEl,dhdx,dhdy);           % compute kinematic matrix

k=k+kinmtx2'*matmtrx*kinmtx2*wtx*wty*detjacob;      % compute element matrix

end
end                                 % end of numerical integration loop

index=FemEldof(nd,No_nodeEl,No_dof);  % extract system dofs associated with element

kk=FemAsmbl1(kk,k,index);             % assemble element matrices 

end
%----------------------------
%  (8) apply boundary conditions
%----------------------------
[kk,ff]=FemAplyc2(kk,ff,bcdof,bcval);
%----------------------------
%  (9) solve the matrix equation
%----------------------------
displmt=kk\ff;

num=1:1:Sys_dof;
displacement=[num' displmt]           % print nodal displacements
%-------------------------------
%  (10) element stress computation
%-------------------------------
for ielp=1:No_element               % loop for the total number of elements

for i=1:No_nodeEl
nd(i)=nodes(ielp,i);          % extract connected node for (iel)-th element
xcoord(i)=gcoord(nd(i),1);            % extract x value of the node
ycoord(i)=gcoord(nd(i),2);            % extract y value of the node
end
%---------------------------
%  (10.1) numerical integration
%---------------------------
intp=0;
for intx=1:nglx
x=point2(intx,1);                     % sampling point in x-axis
wtx=weight2(intx,1);                  % weight in x-axis
for inty=1:ngly
y=point2(inty,2);                     % sampling point in y-axis
wty=weight2(inty,2) ;                 % weight in y-axis
intp=intp+1;

[shape,dhdr,dhds]=FemIsoq4(x,y);      % compute shape functions and
                                      % derivatives at sampling point
jacob2=FemJacobi2(No_nodeEl,dhdr,dhds,xcoord,ycoord);  % compute Jacobian

detjacob=det(jacob2);                 % determinant of Jacobian
invjacob=inv(jacob2);                 % inverse of Jacobian matrix

[dhdx,dhdy]=FemDeriv2(No_nodeEl,dhdr,dhds,invjacob);  % derivatives w.r.t.
                                                 % physical coordinate
kinmtx2=FemKine2D(No_nodeEl,dhdx,dhdy);           % kinematic matrix
 
index=FemEldof(nd,No_nodeEl,No_dof);  % extract system dofs for the element
%---------------------------------------
%  (10.2) extract element displacement vector
%---------------------------------------
for i=1:El_dof
eldisp(i)=displmt(index(i));
end

kinmtx2=FemKine2D(No_nodeEl,dhdx,dhdy);    % compute kinematic matrix

estrain=kinmtx2*eldisp;               % compute strains
estress=matmtrx*estrain;              % compute stresses
 
for i=1:3
strain(intp,i)=estrain(i);                % store for each element
stress(intp,i)=estress(i);                % store for each element
end

location=[ielp,intx,inty]                % print location for stress
stress(intp,:)                         % print stress values

end
end                                 % end of integration loop

end
%--------------------------------------------------------------------------
%    The end
%--------------------------------------------------------------------------
