%--------------------------------------------------------------------------
function [k,m]=BeamElement11(prop,leng,Opt_mass)
%--------------------------------------------------------------------------
%  Purpose: 
%     To calculate stiffness and mass matrices for the beam element. 
%     (Euler-Bernoulli beam) 
%     nodal dof: {v1  1  v2  2} 
%  Synopsis: 
%     [k,m]=BeamElement11(prop,leng,Opt_mass) 
%  Variable Description: 
%     k, m - element stiffness matrix and element mass matrix 
%     prop - the properties of materials and geometry 
%     leng - element length 
%     Opt_mass  = 1 - consistent mass matrix 
%               = 2 - lumped mass matrix 
%--------------------------------------------------------------------------
%  (0) evaluation of the constants
%--------------------------------------------------------------------------
    E=prop(1);                                               % elastic modulus
    u=prop(2);                                                % Poisson's ratio
    rho=prop(3);                           %  mass density (mass per unit volume)

    A=prop(6);                                      % area of beam cross-section
    Iz=prop(8);                  % 2nd moment of inertia of cross-section about axis z

    G=E/(2*(1+u));                                            % shear modulus
%--------------------------------------------------------------------------
%  (1) element stiffness matrix
%--------------------------------------------------------------------------
    c=E*Iz/(leng^3);
    k0=[12        6*leng     -12        6*leng;
        6*leng    4*leng^2   -6*leng    2*leng^2;
       -12       -6*leng      12       -6*leng;
        6*leng    2*leng^2   -6*leng    4*leng^2];
    k=c*k0;
%--------------------------------------------------------------------------
%  (2) element mass matrix
%--------------------------------------------------------------------------
if Opt_mass==1
%---------------------------
%  (2.1) consistent mass matrix
%---------------------------
   mass=rho*A*leng;
   m0=[156        22*leng     54        -13*leng;
       22*leng    4*leng^2    13*leng   -3*leng^2;
       54         13*leng     156       -22*leng;
      -13*leng   -3*leng^2   -22*leng    4*leng^2];
   m=mass/420*m0;
elseif Opt_mass==2
%--------------------------
%  (2.2) lumped mass matrix
%--------------------------
    mass=rho*A*leng;
    m0=diag([1  0  1  0]);
    m=mass/2*m0;
end
%--------------------------------------------------------------------------
%     The end
%--------------------------------------------------------------------------
