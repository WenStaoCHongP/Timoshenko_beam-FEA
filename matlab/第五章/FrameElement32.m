%--------------------------------------------------------------------------
function [k,m]=FrameElement32(prop,leng,xi,al,Opt_section,Opt_mass)
%--------------------------------------------------------------------------
%  Purpose:
%     To calculate stiffness and mass matrices for the 3-d frame element
%     (Timoshenko beam)
%     nodal dof: {u1  v1  w1   x1   y1   z1  u2  v2  w2   x2   y2   z2}
%  Synopsis:
%     [k,m]=FrameElement32(prop,leng,xi,al,Opt_section,Opt_mass)
%  Variable Description:
%     k, m - element stiffness matrix and element mass matrix
%     prop - the properties of materials and geometry: 
%     leng - element length
%     xi - the first row of the coordinate transform matrix between the local and global axes
%        xi(1)=cos(x,x'), xi(2)=cos(x,y'), xi(3)=cos(x,z')
%     al - angle between the reference coordinate system and the local coordinate system
%        for the space element
%     Opt_section - option for type of cross-section
%        = 1 - rectangular cross-section
%        = 2 - circular cross-section
%     Opt_mass - option for mass matrix
%              = 1 - consistent mass matrix
%              = 2 - lumped mass matrix
%--------------------------------------------------------------------------
% (0) calculation of the constants
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
  d=prop(5);                                % inner diameter of beam cross-section
end
A=prop(6);                                           % area of beam cross-section
Iy=prop(7);                       % 2nd moment of inertia of cross-section about axis y
Iz=prop(8);                       % 2nd moment of inertia of cross-section about axis z

polrmoment=prop(9);                        % selection of the polar moment of inertia
shmodule=prop(10);                                   % selection of shear modulus
%------------------------------------
%  (0.2) calculate polar moment of inertia
%------------------------------------
                      % coefficients of the rectangle torsion bar for torsion deformation
 rt=[1.0, 1.2, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 50.0];
 bt=[0.141, 0.166, 0.196, 0.229, 0.249, 0.263, 0.281, 0.291, 0.299, 0.307, 0.313, 0.333];

if Opt_section==1

  r0=h/b;
  if r0<1.0
    r0=1/r0; hx=h; h=b; b=hx;
  end

  ni=length(rt);
  b0=bt(1);

  if r0>=rt(end)
    b0=0.333;
  else
    for ii=1:ni-1                         % linear interpolation for torsion coefficient
      if r0>=rt(ii)
        b0=bt(ii)+(bt(ii+1)-bt(ii))*(r0-rt(ii))/(rt(ii+1)-rt(ii));
      end
    end
  end

end


if polrmoment==0
  Jx=0;
elseif polrmoment==1
  if Opt_section==1; Jx=b0*h*b^3;
  elseif Opt_section==2; Jx=pi*(D^4-d^4)/32; end
else
  Jx=prop(9);
end
%-------------------------------
%  (0.3) calculate the shear modulus
%-------------------------------
 if shmodule==0
   G=0;
 elseif shmodule==1
   G=E/(2*(1+u));
 else
   G=prop(10);
 end
%----------------------------------------------
%  (0.4) selection of correction factor for shear energy
%----------------------------------------------
if Opt_section==1
%--------------------------------------------------------------------------
   ck=6/5;                         % correction factor of the rectangular cross-section
elseif Opt_section==2
   ck=10/9;                           % correction factor of the circular cross-section
end
%--------------------------------------------------------------------------
%  (1) coordinate system transform matrix
%--------------------------------------------------------------------------
 cc=cos(al); ss=sin(al);
 c1=xi(1);   c2=xi(2);   c3=xi(3);
 la=sqrt(c1^2+c2^2);
 
if la==0
  tt=[ 0,  0,  1; -ss,  cc,  0; -cc, -ss,  0];
else
  d1=-c2/la; d2=c1/la; d3=0;
  e1=-c1*c3/la; e2=-c2*c3/la; e3=la;
 
  t1=[ 1,   0,   0; 0,  cc,  ss; 0, -ss,  cc]; 
  t2=[c1, c2, c3; d1, d2, d3; e1, e2, e3];
 
  tt=t1*t2;
end
 
 t0=zeros(3,3);
 T=[tt t0 t0 t0;...
    t0 tt t0 t0;...
    t0 t0 tt t0;...
    t0 t0 t0 tt];
%--------------------------------------------------------------------------
%  (2) element matrix in local coordinate system
%--------------------------------------------------------------------------
%----------------------
%  (2.1) stiffness matrix
%----------------------
 ka=E*A/leng;   kb=G*Jx/leng;
 kc=E*Iz/leng;  kd=E*Iy/leng;
 ke=G*A/(4*ck*leng);

 k11=[ ka            0           0           0             0             0;
        0        4*ke           0           0             0      2*ke*leng;
        0           0        4*ke           0     -2*ke*leng             0;
        0           0           0          kb             0             0;
        0           0   -2*ke*leng           0   ke*leng^2+kd             0;
        0    2*ke*leng           0           0             0  ke*leng^2+kc];
 k12=[-ka            0           0           0             0             0;
        0        -4*ke           0           0             0      2*ke*leng;
        0           0        -4*ke           0     -2*ke*leng             0;
        0           0           0          -kb             0             0;
        0           0    2*ke*leng           0    ke*leng^2-kd             0;
        0   -2*ke*leng           0           0              0   ke*leng^2-kc];
 k21=k12';
 k22=[ ka            0           0           0             0              0;
        0        4*ke           0           0              0     -2*ke*leng;
        0           0        4*ke           0      2*ke*leng              0;
        0           0           0          kb              0             0;
        0           0    2*ke*leng           0   ke*leng^2+kd              0;
        0   -2*ke*leng           0           0              0   ke*leng^2+kc];
 k0=[k11, k12; k21, k22];

 if Jx==0
   k0(:,3)=0; k0(:,5)=0; k0(:,9)=0; k0(:,11)=0;
 end
%--------------------------------------------
%  (2.2) element matrix in global coordinate system
%--------------------------------------------
 k=T'*k0*T;
%--------------------------------------------------------------------------
%  (3) mass matrix
%--------------------------------------------------------------------------

if Opt_mass==1
%-----------------------------
%  (3.1)  consistent mass matrix
%-----------------------------
  ma=rho*A*leng/6; mb=rho*Jx*leng/6;

  m11=[ 2           0           0           0           0           0;
        0           2           0           0           0           0;
        0           0           2           0           0           0;
        0           0           0     2*mb/ma           0           0;
        0           0           0           0           0           0;
        0           0           0           0           0           0];
  m12=[ 1           0           0           0           0           0;
        0           1           0           0           0           0;
        0           0           1           0           0           0;
        0           0           0       mb/ma           0           0;
        0           0           0           0           0           0;
        0           0           0           0           0           0];
  m21=m12';
  m22=[ 2           0           0           0           0           0;
        0           2           0           0           0           0;
        0           0           2           0           0           0;
        0           0           0     2*mb/ma           0           0;
        0           0           0           0           0           0;
        0           0           0           0           0           0];
  m0=ma*[m11, m12; m21, m22];

else
%-------------------------
%  (3.2) lumped mass matrix
%-------------------------
  mass=rho*A*leng/2;
  mb=rho*Jx*leng/2;
  m0=mass*diag([1, 1, 1, mb/mass, 0, 0,... 
               1, 1, 1, mb/mass, 0, 0]);

end
%--------------------------------------------
%  (3.3) element matrix in global coordinate system
%--------------------------------------------
 m=T'*m0*T;
%--------------------------------------------------------------------------
%     The end
%--------------------------------------------------------------------------
