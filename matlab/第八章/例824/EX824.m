%-------------------------------------------------------------------------------------
% Example 8.2.4  
% optimal control design in the coupling(independent)modal space for space 36-bar piezoelectric 
% intelligent truss                                                                                                                                     
% Problem description                                                        
% the modal control of a 36-bar space intelligent truss as shown in Fig. 8.6           
%                                                                            
% Variable descriptions                                                      
%   k = element stiffness matrix                                             
%   m = element mass matrix
%   kk = system stiffness matrix                                             
%   mm = system mass vector                                                 
%   index = a vector containing system dofs associated with each element     
%   gcoord = global coordinate matrix
%   nodes = nodal connectivity matrix for each element
%   bcdof = a vector containing dofs associated with boundary conditions     
%   bcval = a vector containing boundary condition values associated with    
%           the dofs in 'bcdof'         
%-----------------------------------------------------------------------------------
D =[7    24  26];                 % the actuator positions
R = 1.0e-003 *[0.7987    0.8875    0.7917];  % the control gains

%------------------------------------------
% initialize structural parameters
%-----------------------------------------
nel=36;          % number of elements
nnel=2;          % number of nodes per element
ndof=3;          % number of dofs per node
nnode=12;        % total number of nodes in system
sdof=nnode*ndof;  % total system dofs  

%---------------------------
%  nodal coordinates
%-------------------------
gcoord=[0 0 0;1 0 0;1 1 0;0 1 0;...
      0 0 1;1 0 1;1 1 1;0 1 1;...
      0 0 2;1 0 2;1 1 2;0 1 2];

%-----------------------------
% nodal connectivity
%-----------------------------
nodes=zeros(36,2);
nodes(1,1)=1;    nodes(1,2)=5;
nodes(2,1)=2;    nodes(2,2)=6;
nodes(3,1)=3;    nodes(3,2)=7;
nodes(4,1)=4;    nodes(4,2)=8;
nodes(5,1)=5;    nodes(5,2)=6;
nodes(6,1)=6;    nodes(6,2)=7;
nodes(7,1)=7;    nodes(7,2)=8;
nodes(8,1)=5;    nodes(8,2)=8;
nodes(9,1)=6;    nodes(9,2)=8;
nodes(10,1)=5;   nodes(10,2)=7;
nodes(11,1)=1;   nodes(11,2)=8;
nodes(12,1)=4;   nodes(12,2)=5;
nodes(13,1)=4;   nodes(13,2)=7;
nodes(14,1)=3;   nodes(14,2)=8;
nodes(15,1)=3;   nodes(15,2)=6;
nodes(16,1)=2;   nodes(16,2)=7;
nodes(17,1)=2;   nodes(17,2)=5;
nodes(18,1)=1;   nodes(18,2)=6;

for i=1:1
   for j=1:18
      nodes(i*18+j,1)=nodes(j,1)+4*i;
      nodes(i*18+j,2)=nodes(j,2)+4*i;
   end
end

%---------------------------
%  applied constraints
%---------------------------
for i=1:12
   bcdof(i)=i;
   bcval(i)=0;
end
ncon=length(bcdof);            % the number of constraint degrees of freedom
%----------------------------
%  initialization to zero
%----------------------------
leng=zeros(nel,1);             % element length
beta=zeros(nel,1);             % angle between local and global axes
ff=zeros(sdof,1);              % system force vector for initial disturbance
kk=zeros(sdof,sdof);           % system stiffness matrix
mm=zeros(sdof,sdof);          % system mass matrix 
index=zeros(nnel*ndof,1);       % index vector
k=zeros(nnel*ndof,nnel*ndof);   % element stiffness matrix
m=zeros(nnel*ndof,nnel*ndof);  % element mass matrix
B0=zeros(sdof,nel);            % total possible assigned matrix for piezoelectric actuators
B=zeros((sdof-ncon),nact);      % the actuator distributing matrix 
%----------------------------
%  applied nodal force
%----------------------------

ff(35,1)=10000;          %  12th node has 10000 N in y- direction
%----------------------------
%  loop for elements
%----------------------------
for iel=1:nel    % loop for the total number of elements
    nd(1)=nodes(iel,1);   % 1st connected node for the (iel)-th element
    nd(2)=nodes(iel,2);   % 2nd connected node for the (iel)-th element

    x1=gcoord(nd(1),1); y1=gcoord(nd(1),2); z1=gcoord(nd(1),3); % coordinate of 1st node
    x2=gcoord(nd(2),1); y2=gcoord(nd(2),2); z2=gcoord(nd(2),3); % coordinate of 2nd node

    leng(iel)=sqrt((x2-x1)^2+(y2-y1)^2+(z2-z1)^2);  % the length for (iel)-th element

    c(iel)=sqrt((x2-x1)^2)/leng(iel);
    s(iel)=sqrt((y2-y1)^2)/leng(iel);
    r(iel)=sqrt((z2-z1)^2)/leng(iel);
  
   for i=1:length(D)
      if D(i)==iel
             
     [k,m]=fetruss31(leng(iel),c(iel),s(iel),r(iel));      % compute element stiffness and mass matrix of piezoelectric bar
       else
         el=7.2e10;                           %  elastic modulus
         area=0.0025;                         %  cross-sectional area
         rho=2710;                           %   mass density 
        [k,m]=fetruss3(el,leng(iel),area,rho,c(iel),s(iel),r(iel));
   
      end
   end
   
index=feeldof(nd,nnel,ndof);  % extract system dofs for the element

kk=feasmbl1(kk,k,index);      % assemble into system stiffness matrix
mm=feasmbl1(mm,m,index);      % assemble into system mass matrix

end
%-------------------------------------------------------------------
%  apply constraints and solve the matrix for initial disturbance
%-------------------------------------------------------------------

[kk,ff]=feaplyc2(kk,ff,bcdof,bcval);    % apply the boundary conditions
disp=kk\ff;                        % solve the matrix equation to find nodal displacements
disp=disp(13:sdof);
%-------------------------------------------------------------------
% apply optimal control in modal space 
%-------------------------------------------------------------------
kk1=kk(13:sdof,13:sdof);
mm1=mm(13:sdof,13:sdof);
[V,D1]=eig(kk1,mm1);                            % compute eigenvalue and eigenvector 
[lambda,k]=sort(diag(real((D1))));
V=real(V(:,k));
Factor=diag(V'*mm1*V);
Vnorm=V*inv(sqrt(diag(Factor)));                   % normalize eigenvector
natfreq=sqrt(lambda)/(2*pi);
nt=inv(Vnorm)*disp;
nt=real(nt);

%------------------------------------------------------------------
% determine the actuator distributing matrix,including two step:
%  first step: to determine system distributing matrix
%  second step: extract actuator distributing matrix 
%-----------------------------------------------------------------
% first step
%---------------
    for i=1:nact
       if rem(leng(D(i)),sqrt(2))==1
          B0(3*nodes(D(i),1)-2,D(i))=-0.3875588*c(D(i));
          B0(3*nodes(D(i),1)-1,D(i))=-0.3875588*s(D(i));
          B0(3*nodes(D(i),1),D(i))=-0.3875588*r(D(i));
          B0(3*nodes(D(i),2)-2,D(i))=-0.3875588*c(D(i));
          B0(3*nodes(D(i),2)-1,D(i))=-0.3875588*s(D(i));
          B0(3*nodes(D(i),2),D(i))=-0.3875588*r(D(i));
          
       else
          B0(3*nodes(D(i),1)-2,D(i))=-0.6454421*c(D(i));
          B0(3*nodes(D(i),1)-1,D(i))=-0.6454421*s(D(i));
          B0(3*nodes(D(i),1),D(i))=-0.6454421*r(D(i));
          B0(3*nodes(D(i),2)-2,D(i))=-0.6454421*c(D(i));
          B0(3*nodes(D(i),2)-1,D(i))=-0.6454421*s(D(i));
          B0(3*nodes(D(i),2),D(i))=-0.6454421*r(D(i));

       end
    end
 %---------------
 % second step
 %---------------

for i=1:nact
for j=1:nel
   if D(i)==j
      for k=1:(sdof-ncon)
         B(k,i)=B0((k+ncon),j);
      end
   end
end
end
%----------------------------------------------------------------------------------
% compute modal control force and nodal displacements
%  based on couple mode control
%---------------------------------------------------------------------------------
ncm=6;                                        %  the number of controlled modes
A=[zeros(ncm)  eye(ncm)
   -diag(lambda(1:ncm)) zeros(ncm)];
B1=[zeros(ncm,nact)
   Vnorm(:,1:ncm)'*B];
%-------------------------------------------------------------------------------------
% judge whether the system is controlled
%-------------------------------------------------------------------------------------
indexcon=ctrb(A, B1);
rankc=rank(indexcon); 

   t=0:0.001:1;
   n=1/0.001+1;
   x=zeros((sdof-ncon),n);
   x1=zeros((sdof-ncon),n);
   C=eye(2*ncm);
   D2=0;
 %------------------------------------------------------------------------
 % the initial values of the modal displacements and velocities
 %----------------------------------------------------------------------

   for i=1:ncm
      q(i,1)=nt(i);                         
   end
   for i=(ncm+1):(2*ncm)
      q(i,1)=0.0;
   end
   
       yb=initial(ss(A,B1,C,D2),q,t);      %  the modal state vector before control
   for i=1:ncm
       Q(i,1)=real((lambda(i)));
   end
   for i=(ncm+1):(2*ncm)
       Q(i,1)=1;
   end
       Q=diag(Q);
   R1=diag(R);
   [kc,P]=lqr(A,B1,Q,R1);
   A1=A-B1*kc;   
   yc=initial(ss(A1,B1,C,D2),q,t);
   u=-inv(R1)*B1'*P*yc';                                  % applied volatges 

   %-------------------------------------------------------------------------------------
   %  compute displacement before and control only considering controlled modes
   %-------------------------------------------------------------------------------------
   
   x=Vnorm(:,1:ncm)*yb(:,1:ncm)';            % compute nodal displacement before control
   x1=Vnorm(:,1:ncm)*yc(:,1:ncm)';           % compute nodal displacement after control
   %
%-------------------------------------------------------------------------------------    plot(t,x(23,:),'-.',t, x1(23,:))
legend('before control','after control')
xlabel('Tim(sec)','FontWeight','bold')
ylabel('Displacement(m)','FontWeight','bold')
title('节点12的y方向位移')
figure
plot(t,u(1,:))
xlabel('Tim(sec)','FontWeight','bold')
ylabel('Actuation Voltage(V)','FontWeight','bold')
title('杆号7的作动器')
figure
plot(t,u(2,:))
xlabel('Tim(sec)','FontWeight','bold')
ylabel('Actuation Voltage(V)','FontWeight','bold')
title('杆号24的作动器')
figure
plot(t,u(3,:))
xlabel('Tim(sec)','FontWeight','bold')
ylabel('Actuation Voltage(V)','FontWeight','bold')
title('杆号26的作动器')
