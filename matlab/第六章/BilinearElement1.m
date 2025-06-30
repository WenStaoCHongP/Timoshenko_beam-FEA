function varargout=BilinearElement1(varargin)
%--------------------------------------------------------------------------
% Purpose:
% This function is used to calculate element stiffness and mass matrixes 
% of bilinear triangular element 
% Synopsis:
%      k=BilinearElement11(Prop, No_nel,No_dof,gcoord,thickness,iopt)
%     [k,m]=BilinearElement1(Prop, No_nel,No_dof,gcoord,thickness,iopt,Opt_mass)
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
if nargin<6 |nargin>7
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
           nglx=2; ngly=2;                   % 2*2 Gauss-Legendre quadrature
           nglxy=nglx*ngly;                  % number of sampling points per element
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
% numerical integratrion for  element stiffness matrix
%--------------------------------------------------------------------------
        point(1)=-0.577350269189626;                % Sampling pointrs & weights
        point(2)=-point(1);
        weight(1)=1.0;
        weight(2)=weight(1);
        for i=1:No_nel
            xcoord(i)=xycoord(i,1);                % extract x value of the node
            ycoord(i)=xycoord(i,2);                % extract y value of the node
        end
        for intx=1:nglx
            x=point(intx);                         % sampling point in x-axis
            wtx=weight(intx);                      %  weight in x-axis
            for inty=1:ngly
                y=point(inty);                     % sampling point in y-axis
                wty=weight(inty);                  % weight in y-axis
                     dhdr(1)=-0.25*(1-y);          % Compute the derivatives of shape functions
                     dhdr(2)=0.25*(1-y);           % at sampling point
                     dhdr(3)=0.25*(1+y);
                     dhdr(4)=-0.25*(1+y);
                     dhds(1)=-0.25*(1-x);
                     dhds(2)=-0.25*(1+x);
                     dhds(3)=0.25*(1+x);
                     dhds(4)=0.25*(1-x);
           
           jacob=zeros(2,2); 
           for i=1:No_nel                          % compute Jacobian
               jacob(1,1)=jacob(1,1)+dhdr(i)*xcoord(i);
               jacob(1,2)=jacob(1,2)+dhdr(i)*ycoord(i);
               jacob(2,1)=jacob(2,1)+dhds(i)*xcoord(i);
               jacob(2,2)=jacob(2,2)+dhds(i)*ycoord(i);
           end
           detjacob=det(jacob);                    % determinant of Jacobian
           invjacob=inv(jacob);                    % inverse of Jacobian matrix
            for i=1:No_nel
               dhdx(i)=invjacob(1,1)*dhdr(i)+invjacob(1,2)*dhds(i);
               dhdy(i)=invjacob(2,1)*dhdr(i)+invjacob(2,2)*dhds(i);
            end
           for i=1:No_nel
               i1=(i-1)*2+1;
               i2=i1+1;
               kinmtx(1,i1)=dhdx(i);
               kinmtx(2,i2)=dhdy(i);
               kinmtx(3,i1)=dhdy(i);
               kinmtx(3,i2)=dhdx(i);
           end

              k=k+kinmtx'*matmtrx*kinmtx*wtx*wty*detjacob*thickness;
            end
        end
 %---------------------------------------
 % output element stiffness matrix
 %---------------------------------------
   varargout{1}=k;
   
       case 7
           Prop=varargin{1};
           No_nel=varargin{2};
           No_dof=varargin{3};
           xycoord=varargin{4};
           thickness=varargin{5};
           iopt=varargin{6};
           Opt_mass=varargin{7};
           nglx=2; ngly=2;                   % 2*2 Gauss-Legendre quadrature
           nglxy=nglx*ngly;                  % number of sampling points per element
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
%--------------------------------------------------------------------------
% numerical integratrion for  element stiffness matrix
%--------------------------------------------------------------------------
  
        point(1)=-0.577350269189626;                % Sampling pointrs & weights
        point(2)=-point(1);
        weight(1)=1.0;
        weight(2)=weight(1);
        for i=1:No_nel
            xcoord(i)=xycoord(i,1);                % extract x value of the node
            ycoord(i)=xycoord(i,2);                % extract y value of the node
        end
        for intx=1:nglx
            x=point(intx);                         % sampling point in x-axis
            wtx=weight(intx);                      %  weight in x-axis
            for inty=1:ngly
                y=point(inty);                     % sampling point in y-axis
                wty=weight(inty);                  % weight in y-axis
                     dhdr(1)=-0.25*(1-y);          % Compute the derivatives of shape functions
                     dhdr(2)=0.25*(1-y);           % at sampling point
                     dhdr(3)=0.25*(1+y);
                     dhdr(4)=-0.25*(1+y);
                     dhds(1)=-0.25*(1-x);
                     dhds(2)=-0.25*(1+x);
                     dhds(3)=0.25*(1+x);
                     dhds(4)=0.25*(1-x);
           
           jacob=zeros(2,2); 
           for i=1:No_nel                          % compute Jacobian
               jacob(1,1)=jacob(1,1)+dhdr(i)*xcoord(i);
               jacob(1,2)=jacob(1,2)+dhdr(i)*ycoord(i);
               jacob(2,1)=jacob(2,1)+dhds(i)*xcoord(i);
               jacob(2,2)=jacob(2,2)+dhds(i)*ycoord(i);
           end
           detjacob=det(jacob);                    % determinant of Jacobian
           invjacob=inv(jacob);                    % inverse of Jacobian matrix
            for i=1:No_nel
               dhdx(i)=invjacob(1,1)*dhdr(i)+invjacob(1,2)*dhds(i);
               dhdy(i)=invjacob(2,1)*dhdr(i)+invjacob(2,2)*dhds(i);
            end
           for i=1:No_nel
               i1=(i-1)*2+1;
               i2=i1+1;
               kinmtx(1,i1)=dhdx(i);
               kinmtx(2,i2)=dhdy(i);
               kinmtx(3,i1)=dhdy(i);
               kinmtx(3,i2)=dhdx(i);
           end

              k=k+kinmtx'*matmtrx*kinmtx*wtx*wty*detjacob;
            end
        end
 %-------------------------------------------------------------------------
 % element mass matrix
 %-------------------------------------------------------------------------
 area=0.5*((xcoord(1)*ycoord(2)+xcoord(2)*ycoord(3)+xcoord(3)*ycoord(1)-  ...
     xcoord(1)*ycoord(3)-xcoord(2)*ycoord(1)-xcoord(3)*ycoord(2))+ ...
     (xcoord(1)*ycoord(3)+xcoord(3)*ycoord(4)+xcoord(4)*ycoord(1)-  ...
     xcoord(1)*ycoord(4)-xcoord(3)*ycoord(1)-xcoord(4)*ycoord(3)));
 if Opt_mass==1                                   % consistent mass matrix
     m=Prop(3)*thickness*area/36*[4 0 2 0 1 0 2 0; ...
         0 4 0 2 0 1 0 2; ...
         2 0 4 0 2 0 1 0;  ;...
         0 2 0 4 0 2 0 1;  ...
         1 0 2 0 4 0 2 0;  ...
         0 1 0 2 0 4 0 2;  ...
         2 0 1 0 2 0 4 0;  ...
         0 2 0 1 0 2 0 4];
 else                                             % lumped mass matrix
     m=Prop(3)*thickness*area/4*eye(8);
 end
 %-------------------------------------------------------------------------------------
 % output element stiffness matrix and element mass matrix
 %-------------------------------------------------------------------------------------
     varargout{1}=k;
     varargout{2}=m;
   end
end  
%----------------------------------------------------------
%  The end
%----------------------------------------------------------
