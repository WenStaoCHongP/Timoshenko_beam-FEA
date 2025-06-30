%--------------------------------------------------------------------------
function [eta,y,omega1,sdof2]=HarmonicRespt(kk,mm,fd,omega0,t,C,q0,dq0,a,b)
%--------------------------------------------------------------------------
%  Purpose:
%    The function subroutine HarmonicRespt.m calculates harmonic response for a damping
%    or undamping structural system using modal analysis. It uses modal coordinate equations
%    to compute modal responses analytically, then convert the modal responses into
%    physical responses through the coordinate transformation.
%  Synopsis:
%    [eta,y,omega1,sdof2]= HarmonicRespt (kk,mm,fd,omega0,t,C,q0,dq0,a,b)
%  Variable Description:
%    Input parameters:
%      kk, mm - System stiffness and mass matrices
%      fd - Input or forcing influence matrix
%      omega0 - Frequency of the excitation
%      t - Time of evaluation
%      C - Output matrix
%      q0, dq0 - Initial conditions
%      a, b - Parameters for proportional damping [C]=a[M]+b[K]
%    Output parameters:
%      eta - modal coordinate response
%      y - physical coordinate response
%      omega1 - natural frequency
%--------------------------------------------------------------------------
%  (1) Solve the eigenvalue problem
%--------------------------------------------------------------------------
t=t';
[sdof,n1]=size(kk);
[nstep,n2]=size(t);

[V,D]=eig(kk,mm);                       % compute the eigenvalues and eigenvectors
[lambda,ki]=sort(diag(D));                     % sort the eigenvalues and eigenvectors
omega=sqrt(lambda);                                        % natural frequencies
omega1=sqrt(lambda)/(2*pi);                            % the frequency vector in Hz

V=V(:,ki); 
%--------------------------------------------------------------------------
%  (2) Check the eigenvalues
%--------------------------------------------------------------------------
                            % check whether the eigenvalues are infinite and eliminate
                               % the eigenvectors associated with the bad eigenvalues
jk=0;

for ii=1:sdof                                    % loop for find the infinite in omega

  check=omega(ii);
  if check>1.0e12
    jk=jk+1;                                    % location of the infinite frequency
    omi(jk)=ii;                         % storing the location of the infinite frequency
  end

end

sdof2=sdof-jk;
V1=[V(:,1:sdof2)];                                     % truncate the modal vectors
%--------------------------------------------------------------------------
%  (3) Normalize the eigenvectors and compute parameters
%--------------------------------------------------------------------------
Factor=diag(V1'*mm*V1);
Vnorm=V1*inv(sqrt(diag(Factor)));                    %  eigenvectors are normalized
omega2=diag(sqrt(Vnorm'*kk*Vnorm));                         % natural frequencies
           
Fnorm=Vnorm'*fd;                                      % modal input force vector
%--------------------------------------------------------------------------
%  (2) Compute modal damping matrix from the proportional damping matrix
%--------------------------------------------------------------------------
Modamp=Vnorm'*(a*mm+b*kk)*Vnorm;                % form the Rayleigh damping
zeta=diag((1/2)*Modamp*inv(diag(omega2)));                    % the damping ratio

if (max(zeta) >= 1),
  disp('Warning - Your maximum damping ratio is grater than or equal to 1')
  disp('You have to reselect a and b ')
  pause
  disp('If you want to continue, type return key')
end
%--------------------------------------------------------------------------
%  (3) Find out harmonic response of each modal coordinate analytically
%--------------------------------------------------------------------------
eta0=Vnorm'*mm*q0; deta0=Vnorm'*mm*dq0;
                % initial conditions for modal coordinates both displacement and velocity
eta=zeros(nstep,sdof2);

phase0=omega0*t;

for i=1:sdof2                                  % responses are obtained for n modes
  gama=omega0/omega(i);
  omegad=omega(i)*sqrt(1-zeta(i)^2);
  phase=omegad*t;
  Exx=exp(-zeta(i)*omega(i)*t);

  C1=eta0(i);
  C2=(deta0(i)+eta0(i)*zeta(i)*omega(i))/omegad;

  X0=sqrt((1-gama^2)^2+(2*zeta(i)*gama)^2);
  XX=Fnorm(i)/(omega(i)^2*X0);
  XP=atan((2*zeta(i)*gama)/(1-gama^2));

  D1=(zeta(i)*omega(i)*cos(XP)+omega0*sin(XP))/omegad;
  D2=cos(XP);

  eta(:,i)=C1*Exx.*cos(phase)+C2*Exx.*sin(phase)...
           -XX*Exx.*(D1*sin(phase)+D2*cos(phase))...
           +XX*cos(phase0-XP);
                         % response included transient-state for the harmonic excitation
end
%--------------------------------------------------------------------------
%  (4) Convert modal coordinate responses to physical coordinate responses
%--------------------------------------------------------------------------
eta=eta';
y=C*Vnorm*eta;
 
if (a+b)==0
  disp('The response results of undamping system')
else
  disp('The response results of damping system')
end

disp('The excitation is harmonic force')
%--------------------------------------------------------------------------
%     The end
%--------------------------------------------------------------------------
