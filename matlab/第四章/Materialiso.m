%--------------------------------------------------------------------------
function [matmtrx]=Materialiso(iopt,E,u)
%--------------------------------------------------------------------------
%  Purpose:
%     determine the constitutive equation for isotropic material
%  Variable Description:
%     E - elastic modulus
%     u - Poisson's ratio
%     iopt=1 - plane stress analysis
%     iopt=2 - plane strain analysis
%     iopt=3 - axisymmetric analysis
%     iopt=4 - three dimensional analysis
%--------------------------------------------------------------------------
 if iopt==1        % plane stress
   matmtrx= E/(1-u*u)* ...
   [1  u   0; ...
   u   1   0; ...
   0   0  (1-u)/2];
 elseif iopt==2     % plane strain
   matmtrx= E/((1+u)*(1-2*u))* ...
   [(1-u)  u      0; 
   u    (1-u)     0;
   0      0   (1-2*u)/2];
 elseif iopt==3     % axisymmetry
   matmtrx= E/((1+u)*(1-2*u))* ...
   [(1-u)  u     u     0; 
   u    (1-u)    u     0;
   u      u   (1-u)    0;
   0      0     0   (1-2*u)/2];
 else     % three-dimension
   matmtrx= E/((1+u)*(1-2*u))* ...
   [(1-u)  u     u     0       0       0; 
   u     (1-u)   u     0       0       0;
   u      u    (1-u)   0       0       0;
   0      0     0    (1-2*u)/2  0       0;
   0      0     0     0      (1-2*u)/2  0;
   0      0     0     0       0       (1-2*u)/2];
 end
%--------------------------------------------------------------------------
%    The end
%--------------------------------------------------------------------------
