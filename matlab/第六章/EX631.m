%------------------------------------------------------------
% Example 6.3.1  
% Dynamic characterisitic analysis and static analysis of the
% plate for plane stress analysis.
%
% Variable descriptions
%   dk = element stiffness matrix
%   dm = element mass matrix
%    k = system stiffness matrix 
%    m = system mass matrix 
%    gcoord = coordinate values of each node
%    nodes = nodal connectivity of each element
%гнгнгнгнгнгнгнгнгнгнгнгнгнгнгнгнгнгнгнгнгнгнгнгнгнгнгнгнгнгнгнгн
% Input data for control paramters
%-----------------------------------------------------------

Element_number=100;       % number of elements
No_nel=3;                 % number of nodes per element
No_dof=2;                 % number of dofs per node
Node_number=66;          % total number of nodes in system
Prop(1)=6.9e10;           % elastic modulus
Prop(2)=0.3;              % Poisson's ratio
Prop(3)=2700;             % Mass density
t=0.001;                  % Thickness
%----------------------------------------------------------
% Input data for nodal coordinate values
% gcoord(i,j) where i-> node no. and j-> x or y
%----------------------------------------------------------
for loopi=1:Node_number
    if rem(loopi,11)~=0
        gcoord(loopi,1)=(loopi-floor(loopi/11)*11-1)*0.02;
        gcoord(loopi,2)=floor(loopi/11)*0.02;
    else
        gcoord(loopi,1)=0.2;
        gcoord(loopi,2)=(floor(loopi/11)-1)*0.02;
    end
end
%----------------------------------------------------------
% Input data for nodal connectivity for each element 
% nodes(i,j) where i-> element no. and j-> connected nodes
%--------------------------------------------------------
No_element=0;
for loopi=1:10
    for loopj=1:10
        No_element=No_element+1;
       if rem(No_element,2)~=0
          nodes(No element,1)=(No_element+1)/2+floor((loopi-1)/2);
          nodes(No_element,2)=(No_element+1)/2+1+floor((loopi-1)/2);
          nodes(No_element,3)=(No_element+1)/2+11+floor((loopi-1)/2); 
       else
          nodes(No_element,1)=(No_element+2)/2+floor((loopi-1)/2);
          nodes(No_element,2)=(No_element+2)/2+11+floor((loopi-1)/2);
          nodes(No_element,3)=(No_element+2)/2+10+floor((loopi-1)/2); 
       end
    end
end 
%------------------------------------------------------
% Inout data for boundary conditions
%------------------------------------------------------
ed(1:Node_number,1:2)=1;            %elment_displacement 
constraint=[1,1;1,2;12,1;12,2;23,1;23,2;34,1;34,2;45,1;45,2;56,1;56,2];
 
for loopi=1:length(constraint);
    ed(constraint(loopi,1),constraint(loopi,2))=0;
end
dof=0;
for loopi=1:Node_number
     for loopj=1:2
        if ed(loopi,loopj)~=0
            dof=dof+1;
            ed(loopi,loopj)=dof;
        end
    end
end
%--------------------------------------------------------
%Initilization of matrices 
%--------------------------------------------------------
k=zeros(dof,dof);          % system stiffness matrix
m=zeros(dof,dof);          % system mass matrix
f=zeros(Node_number*No_dof,1); % system force vector
disp=zeros(dof,1);         % system displacement vector
eldisp=zeros(No_nel*No_dof,1); % element diaplacement vector
stress=zeros(Element_number,3);  % matrix containing stress components
strain=zeros(Element_number,3);  % matrix containing strain components
kinmtx=zeros(3,No_nel*No_dof);   % kinematic matrix 
matmtrx=zeros(3,3);               % constitutive matrix
e2s(1:6)=0;    % index of transform the elment displament number to structural
%---------------------------
% force vector
%---------------------------
f(21)=5000;
f(131)=5000;
%--------------------------------------------------------
% Compute system stiffness and mass matrices
%-------------------------------------------------------
 for loopi=1:Element_number            %  loop for the total number of elements
     for zi=1:3
        e2s((zi-1)*2+1)=ed(nodes(loopi,zi),1);
        e2s((zi-1)*2+2)=ed(nodes(loopi,zi),2);
     end
   for loopj=1:3
       xycoord(loopj,1)=gcoord(nodes(loopi,loopj),1);
       xycoord(loopj,2)=gcoord(nodes(loopi,loopj),2);
   end
  iopt=1;                               % plane stress analysis
 % iopt=2;                              % plane strain analysis
  Opt_mass=1;                           % Consistent mass matrix
  [dk,dm]=LineartriElement1(Prop,No_nel,No_dof,xycoord,t,iopt,Opt_mass);
    for jx=1:6
        for jy=1:6
            if(e2s(jx)*e2s(jy)~=0)
                k(e2s(jx),e2s(jy))=k(e2s(jx),e2s(jy))+dk(jx,jy);
                m(e2s(jx),e2s(jy))=m(e2s(jx),e2s(jy))+dm(jx,jy);
            end
        end
    end
 end
 %--------------------------------------
 % compute system displacement
 %--------------------------------------
 Number_con=-1;                % apply boundary conditions
 for loopi=1:length(constraint)   
     Number_con= Number_con+1;
     f((constraint(loopi,1)-1)*No_dof+constraint(loopi,2)-Number_con)=[];
 end
 disp=k\f;                    %solve the matrix equation
 %----------------------------------------
 % element stress and strain computation
 %---------------------------------------
 for loopi=1:Element_number 
     for zi=1:3                       % extract system dofs for the element
         e2s((zi-1)*2+1)=ed(nodes(loopi,zi),1);
         e2s((zi-1)*2+2)=ed(nodes(loopi,zi),2);
     end
     for loopj=1:3
        xycoord(loopj,1)=gcoord(nodes(loopi,loopj),1);
        xycoord(loopj,2)=gcoord(nodes(loopi,loopj),2);
     end
     for loopk=1:No_nel*No_dof % extract element displacement vector
         if e2s(loopk)==0
            eldisp(loopk)=0;
         else
            eldisp(loopk)=disp(e2s(loopk));
         end
     end

      x1=xycoord(1,1);y1=xycoord(1,2);
      x2=xycoord(2,1);y2=xycoord(2,2);
      x3=xycoord(3,1);y3=xycoord(3,2);
      area=0.5*(x1*y2+x2*y3+x3*y1-x1*y3-x2*y1-x3*y2);      % area of triangule
      dhdx=(1/(2*area))*[(y2-y3) (y3-y1) (y1-y2)];         % derivatives w.r.t x
      dhdy=(1/(2*area))*[(x3-x2) (x1-x3) (x2-x1)];         % derivatives w.r.t y
      for i=1:No_nel                                       % kinematic matrix
          i1=(i-1)*2+1;
          i2=i1+1;
          kinmtx(1,i1)=dhdx(i);
          kinmtx(2,i2)=dhdy(i);
          kinmtx(3,i1)=dhdy(i);
          kinmtx(3,i2)=dhdx(i);
      end
       matmtrx=Prop(1)/(1-Prop(2)*Prop(2))* ...
                   [1 Prop(2) 0; ...
                   Prop(2) 1 0;  ...
                   0 0 (1-Prop(2))/2];
       estrain=kinmtx*eldisp;                         % compute strains
       estress=matmtrx*estrain ;                      % compute stresses
       for i=1:3
           strain(loopi,i)=estrain(i);               % store for each element
           stress(loopi,i)=estress(i);               % store for each element
       end
   
 end
 %------------------------------------------------------
 % Construct the taotal displacement
 %------------------------------------------------------
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
     if length(displ)>=Node_number*No_dof
             flag=0;
     end 
 end
 %--------------------------------------------
 % Compute the first ten natrual frequencies
 %--------------------------------------------
 p=10;
 epsilon=1e-5;
 [v,d]=Ssiter(k,m,p,epsilon);
 d=sqrt(d)/(2*pi);
 %---------------------------------------------------
 %  print fem solutions
 %---------------------------------------------------
  num=1:1:p;
  frequency=[num' d]                % print frequency
 %
  num=1:1:Node_number*No_dof;
  displace=[num' dispt']            % print nodal displacements
 %
  num=1:Element_number;
stresses=[num' stress]            % print stresses

>> format long

frequency =
  1.0e+004 *
   0.00010000000000   0.18604065704136
   0.00020000000000   0.63693976247305
   0.00030000000000   0.70795736724809
   0.00040000000000   1.49743000911713
   0.00050000000000   1.87501877728104
   0.00060000000000   1.88424652438432
   0.00070000000000   2.42279455024536
   0.00080000000000   2.54380218643713
   0.00090000000000   2.69029118728467
   0.00100000000000   2.71490380740973
displace =
  1.0e+002 *
   0.01000000000000                  0
   0.02000000000000                  0
   0.03000000000000   0.00000034291373
   0.04000000000000   0.00000015963419
    бнбнбн
   1.29000000000000   0.00000327251912
   1.30000000000000  -0.00000032693461
   1.31000000000000   0.00000380946666
   1.32000000000000  -0.00000054632754
stresses =
  1.0e+008 *
   0.00000001000000   1.30005754768989   0.39001726430697   0.21182228466200
   0.00000002000000   1.12293677324908   0.07616439071518   0.07963707173014
   0.00000003000000   1.17720848679054   0.09244590477762   0.04321198916920
   0.00000004000000   1.17304599196842   0.02612618266368   0.02423600768642
   бнбн..
   0.00000097000000   1.17311067119685   0.00959355874785   0.06554658062287
   0.00000098000000   1.52792274071915   0.07363651493596   0.17831524139106
   0.00000099000000   1.14515705809605  -0.04119318985097   0.08300345141375
   0.00000100000000   1.78924143352395  -0.21075856647605   0.21075856647604
