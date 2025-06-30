%--------------------------------------------------------------------------
function [k,m]=FrameElement22(prop,leng,beta,Opt_section,Opt_mass)
%--------------------------------------------------------------------------
%  Purpose:
%     To calculate stiffness and mass matrices for the 2-d frame element
%     (Timoshenko beam)
%     nodal dof: {u1  v1   1  u2  v2   2}
%  Synopsis:
%     [k,m]=FrameElement22(prop,leng,beta,Opt_section,Opt_mass)
%  Variable Description:
%     k, m - element stiffness matrix and element mass matrix
%     prop - the properties of materials and geometry
%     leng - element length
%     beta - angle between the local and global axes
%          is positive if the local axis is in the ccw direction from the global axis
%     Opt_section - option for type of cross-section
%               = 1 - rectangular cross-section
%               = 2 - circular cross-section
%     Opt_mass - option for mass matrix
%             = 1 - consistent mass matrix
%             = 2 - lumped mass matrix
%--------------------------------------------------------------------------
%  (0) calculation of the constants
%--------------------------------------------------------------------------
%-------------------------------
%  (0.1) evaluation of the constants
%-------------------------------
E=prop(1);                                                   % elastic modulus
u=prop(2);                                                    % Poisson's ratio
rho=prop(3);                                %  mass density (mass per unit volume)

if Opt_section==1
  h=prop(4);                                       % height of beam cross-section
  b=prop(5);                                       % width of beam cross-section
elseif Opt_section==2
  D=prop(4);                                % outer diameter of beam cross-section
  d=prop(5);                                 % inner diameter of beam cross-section
end

A=prop(6);                                           % area of beam cross-section
Iy=prop(7);                       % 2nd moment of inertia of cross-section about axis y
Iz=prop(8);                       % 2nd moment of inertia of cross-section about axis z

polrmoment=prop(9);                         % selection of the polar moment of inertia
shmodule=prop(10);                                    % selection of shear modulus
%-------------------------------
%  (0.2) calculate the shear modulus
%-------------------------------
 if shmodule==0
   G=0;
 elseif shmodule==1
   G=E/(2*(1+u));
 else
   G=prop(10);
 end
%----------------------------------------------
%  (0.3) selection of correction factor for shear energy
%----------------------------------------------
if Opt_section==1
   ck=6/5;                         % correction factor of the rectangular cross-section
elseif Opt_section==2
   ck=10/9;                          % correction factor of the circular cross-section
end 
%--------------------------------------------------------------------------
%  (1) rotation matrix for the coordinate transformation
%--------------------------------------------------------------------------
cc=cos(beta); ss=sin(beta);

 T=[ cc   ss    0    0    0    0;
    -ss   cc    0    0    0    0;
     0    0    1    0    0    0;
     0    0    0    cc   ss    0;
     0    0    0   -ss   cc    0;
     0    0    0    0    0    1];
%--------------------------------------------------------------------------
%  (2) stiffness matrix
%--------------------------------------------------------------------------
%----------------------------------
%  (2.1) stiffness matrix at the local axis
%----------------------------------
 ka=E*A/leng;
 kc=E*Iz/leng;  
 kd=G*A/(4*ck*leng);
 
 k0=[ka    0          0            -ka     0           0;
      0   4*kd        2*kd*leng      0    -4*kd        2*kd*leng;
      0   2*kd*leng    kc+kd*leng^2   0    -2*kd*leng   -kc+kd*leng^2;
    -ka    0           0             ka    0           0;
      0  -4*kd        -2*kd*leng      0     4*kd       -2*kd*leng;
      0   2*kd*leng   -kc+kd*leng^2   0     -2*kd*leng   kc+kd*leng^2];
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
    mass=rho*A*leng;
    m0=mass/6*[2     0     0     1     0     0;
               0     2     0     0     1     0;
               0     0     0     0     0     0;
               1     0     0     2     0     0;
               0     1     0     0     2     0;
               0     0     0     0     0     0];
else
%-------------------------
%  (3.2) lumped mass matrix
%-------------------------
    mass=rho*A*leng;
    m0=mass/2*diag([1, 1, 0, 1, 1, 0]);
end
%------------------------------
%  (3.3) mass in the global system
%------------------------------
    m=T'*m0*T;
%--------------------------------------------------------------------------
%     The end
% -------------------------------------------------------------------------
