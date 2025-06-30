function varargout=LineartriElement1(varargin)
%--------------------------------------------------------------------------
% Purpose:
% This function is used to calculate element stiffness and mass matrixes 
% of linear triangular element 
% Synopsis:
%      k=LineartriElement1(Prop, No_nel,No_dof,xycoord,thickness,iopt)
%     [k,m]=LineartriElement1(Prop, No_nel,No_dof,xycoord,thickness,iopt,Opt_mass)
% Variable Description:
%       Input parameters
%          Prop(1)- elastic modulus
%          Prop(2)- Poisson's ratio
%          Prop(3)- mass density
%          No_nel - number of nodes per element
%          No_dof - number of dofs per node
%          xycoord -coord values of nodes
%          iopt=1 - plane stress analysis
%          iopt=2 - plane strain analysis
%          thickness - element thickness 
%          Opt_mass =1 - consistent mass matrix 
%          Opt_mass= 2 - lumped mass matrix 
%       Output parameters
%           k - element stiffness matrix
%           m - element mass matrix
% Author: Dr.XU Bin, Time:2006-12-08
%--------------------------------------------------------------------------

if nargin<6 & nargin>7
    error('Incorrect number of input arguments')
else
   switch nargin
       case 6
           Prop=varargin{1};
           No_nel=varargin{2};
           No_dof=varargin{3};
           xycoord=varargin{4};
           thickness=varargin{5};
           iopt=varargin{6};
           k=zeros(No_nel*No_dof);
           kinmtx=zeros(3,No_nel*No_dof);
           matmtrx=zeros(3,3);    
 %-------------------------------------------------------------------------
 % the constitutive equation for isotropic material
 %-------------------------------------------------------------------------
           if iopt==1                                           % constitutive matrix for plane stress
               matmtrx=Prop(1)/(1-Prop(2)*Prop(2))* ...
                   [1 Prop(2) 0; ...
                   Prop(2) 1 0;  ...
                   0 0 (1-Prop(2))/2];
           else
               matmtrx=Prop(1)/((1+Prop(2))*(1-2*Prop(2)))* ...   % constitutive matrix for plane strain
                   [1-Prop(2) Prop(2) 0; ...
                   Prop(2) 1-Prop(2) 0;  ...
                   0 0 (1-2*Prop(2))/2];
           end
%--------------------------------------------------------------------------
% the kinematic equation between strains and displacements 
%--------------------------------------------------------------------------
           x1=xycoord(1,1);y1=xycoord(1,2);
           x2=xycoord(2,1);y2=xycoord(2,2);
           x3=xycoord(3,1);y3=xycoord(3,2);
           area=0.5*(x1*y2+x2*y3+x3*y1-x1*y3-x2*y1-x3*y2);      % area of triangule
           dhdx=(1/(2*area))*[(y2-y3) (y3-y1) (y1-y2)];         % derivatives w.r.t x
           dhdy=(1/(2*area))*[(x3-x2) (x1-x3) (x2-x1)];         % derivatives w.r.t y
           for i=1:No_nel
               i1=(i-1)*2+1;
               i2=i1+1;
               kinmtx(1,i1)=dhdx(i);
               kinmtx(2,i2)=dhdy(i);
               kinmtx(3,i1)=dhdy(i);
               kinmtx(3,i2)=dhdx(i);
           end
 %------------------------------------------------------------------------
 % element stiffness matrix
 %------------------------------------------------------------------------
   k=kinmtx'*matmtrx*kinmtx*area*thickness;

 %---------------------------------------
 % output element stiffness matrix
 %---------------------------------------
   varargout{1}=k;
 %---------------------------------------
   
       case 7
           Prop=varargin{1};
           No_nel=varargin{2};
           No_dof=varargin{3};
           xycoord=varargin{4};
           thickness=varargin{5};
           iopt=varargin{6};
           Opt_mass=varargin{7};
           k=zeros(No_nel*No_dof);
           m=zeros(No_nel*No_dof);
           kinmtx=zeros(3,No_nel*No_dof);
           matmtrx=zeros(3,3);    
 %-------------------------------------------------------------------------
 % the constitutive equation for isotropic material
 %-------------------------------------------------------------------------
           if iopt==1                                           % constitutive matrix for plane stress
               matmtrx=Prop(1)/(1-Prop(2)*Prop(2))* ...
                   [1 Prop(2) 0; ...Prop(2)
                   Prop(2) 1 0;  ...
                   0 0 (1-Prop(2))/2];
           else
               matmtrx=Prop(1)/((1+Prop(2))*(1-2*Prop(2)))* ...   % constitutive matrix for plane strain
                   [1-Prop(2) Prop(2) 0; ...
                   Prop(2) 1-Prop(2) 0;  ...
                   0 0 (1-2*Prop(2))/2];
           end
%----------------------------------------------------------------------------------------
% the kinematic equation between strains and displacements 
%---------------------------------------------------------------------------------------
           x1=xycoord(1,1);y1=xycoord(1,2);
           x2=xycoord(2,1);y2=xycoord(2,2);
           x3=xycoord(3,1);y3=xycoord(3,2);
           area=0.5*(x1*y2+x2*y3+x3*y1-x1*y3-x2*y1-x3*y2);      % area of triangule
           dhdx=(1/(2*area))*[(y2-y3) (y3-y1) (y1-y2)];         % derivatives w.r.t x
           dhdy=(1/(2*area))*[(x3-x2) (x1-x3) (x2-x1)];         % derivatives w.r.t y
           for i=1:No_nel
               i1=(i-1)*2+1;
               i2=i1+1;
               kinmtx(1,i1)=dhdx(i);
               kinmtx(2,i2)=dhdy(i);
               kinmtx(3,i1)=dhdy(i);
               kinmtx(3,i2)=dhdx(i);
           end
 %------------------------------------------------------------------------
 % element stiffness matrix
 %------------------------------------------------------------------------
   k=kinmtx'*matmtrx*kinmtx*area*thickness;   
 
 %-------------------------------------------------------------------------
 % element mass matrix
 %-------------------------------------------------------------------------
 
 if Opt_mass==1                                   % consistent mass matrix
     m=Prop(3)*thickness*area/12*[2 0 1 0 1 0; ...
         0 2 0 1 0 1; ...
         1 0 2 0 1 0; ...
         0 1 0 2 0 1; ...
         1 0 1 0 2 0; ...
         0 1 0 1 0 2];
 else                                             % lumped mass matrix
     m=Prop(3)*thickness*area/3*eye(6);
 end
 %-----------------------------------------------------------
 % output element stiffness matrix and element mass matrix
 %-----------------------------------------------------------
     varargout{1}=k;
     varargout{2}=m;
   end
end
%------------------------------------------------------------
% The end
%------------------------------------------------------------
