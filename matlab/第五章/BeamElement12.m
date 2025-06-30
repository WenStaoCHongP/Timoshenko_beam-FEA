%--------------------------------------------------------------------------
function [k,m]=BeamElement12(prop,leng,Opt_mass)
%--------------------------------------------------------------------------
%  Purpose: 
%     To calculate stiffness and mass matrices for the beam element. 
%     (Timoshenko beam) 
%     nodal dof: {v1  1  v2  2} 
%  Synopsis: 
%     [k,m]=BeamElement12(prop,leng,Opt_mass) 
%  Variable Description: 
%     k, m - element stiffness matrix and element mass matrix 
%     prop - the properties of materials and geometry 
%     leng - element length 
%     Opt_mass  = 1 - consistent mass matrix 
%               = 2 - lumped mass matrix 
%--------------------------------------------------------------------------
%  (0) evaluation of the constants
%--------------------------------------------------------------------------
    E=prop(1);                                                % elastic modulus
    u=prop(2);                                                 % Poisson's ratio
    rho=prop(3);                            %  mass density (mass per unit volume)

    A=prop(6);                                       % area of beam cross-section
    Iz=prop(8);                   % 2nd moment of inertia of cross-section about axis z

    G=E/(2*(1+u));                                             % shear modulus
%--------------------------------------------------------------------------
%  (1) element stiffness matrix
%--------------------------------------------------------------------------
    c=E*Iz/leng;  
    d=(5/6)*G*A/(4*leng);
    k=[4*d         2*d*leng     -4*d         2*d*leng;
       2*d*leng    c+d*leng^2   -2*d*leng   -c+d*leng^2;
      -4*d        -2*d*leng      4*d        -2*d*leng;
       2*d*leng   -c+d*leng^2   -2*d*leng    c+d*leng^2];
%--------------------------------------------------------------------------
%  (2) element mass matrix
%--------------------------------------------------------------------------
if Opt_mass==1
%---------------------------
%  (2.1) consistent mass matrix
%---------------------------
    mass=rho*A*leng;
    m0=[2     0     1     0;
        0     0     0     0;
        1     0     2     0;
        0     0     0     0];
    m=mass/6*m0;
else
%-------------------------
%  (2.2) lumped mass matrix
%-------------------------
    mass=rho*A*leng;
    m0=diag([1  0  1  0]);
    m=mass/2*m0;
end
%--------------------------------------------------------------------------
%     The end
%--------------------------------------------------------------------------
