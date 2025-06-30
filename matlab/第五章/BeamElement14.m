%--------------------------------------------------------------------------
function [k,m]=BeamElement14(prop,leng,Opt_mass)
%--------------------------------------------------------------------------
%  Purpose:
%     To calculate stiffness and mass matrices for the beam element.
%     (mixed beam)
%     nodal dof: {M1  v1  M2  v2}
%  Synopsis:
%     [k,m]=BeamElement14(prop,leng,Opt_mass)
%  Variable Description:
%     k, m - element stiffness matrix and element mass matrix
%     prop - the properties of materials and geometry
%     leng - element length
%     Opt_mass = 1 - consistent mass matrix
%              = 2 - lumped mass matrix
%--------------------------------------------------------------------------
%  (0) evaluation of the constants
%--------------------------------------------------------------------------
    E=prop(1);                                                % elastic modulus
    u=prop(2);                                                 % Poisson's ratio
    rho=prop(3);                            %  mass density (mass per unit volume)
    A=prop(6);                                       % area of beam cross-section
    Iz=prop(8);                   % 2nd moment of inertia of cross-section about axis z
    shmodule=prop(10);                               % selection of shear modulus
%--------------------------------------------------------------------------
%  (1) element stiffness matrix
%--------------------------------------------------------------------------
if shmodule==0
%------------------------------------------------
%  (1.1) stiffness matrix not including shear deformation
%------------------------------------------------
    c=1/(6*E*Iz*leng);
    k=[2*leng^2    6*E*Iz     leng^2      -6*E*Iz;
       6*E*Iz      0         -6*E*Iz      0;
       leng^2     -6*E*Iz     2*leng^2     6*E*Iz;
       -6*E*Iz      0         6*E*Iz      0];
    k=c*k;
elseif shmodule==1
%---------------------------------------------
%  (1.2) stiffness matrix including shear deformation
%---------------------------------------------
    G=E/(2*(1+u));                                            % shear modulus
    c=1/(6*E*Iz*leng);
    d=(6/5)*6*E*Iz/(G*A);
    k=[2*leng^2+d    6*E*Iz     leng^2-d      -6*E*Iz;
       6*E*Iz        0         -6*E*Iz        0;
       leng^2-d     -6*E*Iz      2*leng^2+d    6*E*Iz;
       -6*E*Iz       0          6*E*Iz        0];
    k=c*k;
else
%---------------------------------------------
%  (1.3) stiffness matrix including shear deformation
%---------------------------------------------
    G=prop(10);                                               % shear modulus
    c=1/(6*E*Iz*leng);
    d=(6/5)*6*E*Iz/(G*A);
    k=[2*leng^2+d    6*E*Iz     leng^2-d      -6*E*Iz;
       6*E*Iz        0         -6*E*Iz        0;
       leng^2-d     -6*E*Iz      2*leng^2+d    6*E*Iz;
       -6*E*Iz       0          6*E*Iz        0];
    k=c*k;
end
%--------------------------------------------------------------------------
%  (2) element mass matrix
%--------------------------------------------------------------------------
if Opt_mass==1
%---------------------------
%  (2.1) consistent mass matrix
%---------------------------
    mass=rho*A*leng;
    m0=diag([0  1  0  1]);
    m=mass/2*m0;     
elseif Opt_mass==2
%-------------------------
%  (2.2) lumped mass matrix
%-------------------------
    mass=rho*A*leng;
    m0=diag([0  1  0  1]);
    m=mass/2*m0;
else
%--------------------------
%  (2.3) diagonal mass matrix
%--------------------------
    mass=rho*A*leng;
    m0=diag([1  1  1  1]);
    m=mass/2*m0;
end
%--------------------------------------------------------------------------
%     The end
%--------------------------------------------------------------------------
