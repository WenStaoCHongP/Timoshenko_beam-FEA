%--------------------------------------------------------------------------
function [k,m]=FrameElement21(prop,leng,beta,Opt_section,Opt_mass)
%--------------------------------------------------------------------------
%  Purpose:
%     To calculate stiffness and mass matrices for the 2-d frame element
%     (Euler-Bernoulli beam)
%     nodal dof: {u1  v1   1  u2  v2   2}
%  Synopsis:
%     [k,m]=FrameElement21(prop,leng,beta,Opt_section,Opt_mass)
%  Variable Description:
%     k, m - element stiffness matrix and element mass matrix
%     prop - the properties of materials and geometry:
%     leng - element length
%     beta - angle between the local and global axes
%          is positive if the local axis is in the ccw direction from the global axis
%     Opt_section - option for type of cross-section
%               = 1 - rectangular cross-section
%               = 2 - circular cross-section
%     Opt_mass - option for mass matrix
%             = 1 - consistent mass matrix
%             = 2 - lumped mass matrix
%             = 3 - diagonal mass matrix
%--------------------------------------------------------------------------
%  (0) evaluation of the constants
%--------------------------------------------------------------------------
E=prop(1);                                                    % elastic modulus
u=prop(2);                                                     % Poisson's ratio
rho=prop(3);                                 %  mass density (mass per unit volume)

if Opt_section==1
  h=prop(4);                                        % height of beam cross-section
  b=prop(5);                                        % width of beam cross-section
elseif Opt_section==2
  D=prop(4);                                 % outer diameter of beam cross-section
  d=prop(5);                                 % inner diameter of beam cross-section
end

A=prop(6);                                           % area of beam cross-section
Iy=prop(7);                       % 2nd moment of inertia of cross-section about axis y
Iz=prop(8);                       % 2nd moment of inertia of cross-section about axis z
%--------------------------------------------------------------------------
%  (1) rotation matrix for the coordinate transformation
%--------------------------------------------------------------------------
cc=cos(beta); ss=sin(beta);

 T=[ cc   ss    0    0    0    0;
    -ss   cc    0    0    0    0;
     0    0    1    0    0    0;
     0    0    0   cc   ss     0;
     0    0    0   -ss   cc    0;
     0    0    0    0    0    1];
%--------------------------------------------------------------------------
%  (2) stiffness matrix
%--------------------------------------------------------------------------
%----------------------------------
%  (2.1) stiffness matrix at the local axis
%----------------------------------
 ka=E*A/leng; kc=E*Iz/(leng^3);

 k0=[ka   0            0            -ka     0           0;
      0   12*kc        6*leng*kc      0    -12* kc       6*leng*kc;
      0   6*leng*kc    4*leng^2*kc     0    -6*leng*kc   2*leng^2*kc;
    -ka   0            0             ka     0           0;
      0  -12*kc       -6*leng*kc       0     12*kc      -6*leng*kc;
      0   6*leng*kc    2*leng^2*kc     0    -6*leng*kc    4*leng^2*kc];
%------------------------------------
%  (2.2) stiffness matrix at the global axis
%------------------------------------
 k=T'*k0*T;
%--------------------------------------------------------------------------
%  (3) mass matrix
%--------------------------------------------------------------------------
if Opt_mass==1
%---------------------------
%  (3.1) consistent mass matrix
%---------------------------
  ma=rho*A*leng/6; mb=rho*A*leng/420;

  m0=[2*ma  0           0            ma    0            0;
      0     156*mb      22*leng*mb   0     54*mb       -13*leng*mb;
      0     22*leng*mb  4*leng^2*mb  0     13*leng*mb  -3*leng^2*mb;
      ma    0           0            2*ma  0            0;
      0     54*mb       13*leng*mb   0     156*mb      -22*leng*mb;
      0    -13*leng*mb -3*leng^2*mb  0    -22*leng*mb   4*leng^2*mb];
elseif Opt_mass==2
%-------------------------
%  (3.2) lumped mass matrix
%-------------------------
  mass=rho*A*leng;
  m0=mass/2*diag([1  1  0  1  1  0]);
end
%-----------------------------------
%  (3.3) mass matrix in the global system
%-----------------------------------
m=T'*m0*T;
%--------------------------------------------------------------------------
%     The end
%--------------------------------------------------------------------------
