%--------------------------------------------------------------------------
function [eta,y,omega1,sdof2]=StepRespt(kk,mm,fd,u,t,C,q0,dq0,a,b)
%--------------------------------------------------------------------------
%  Purpose:
%    The function subroutine StepRespt.m calculates step response for a damping or
%    undamping structural system using modal analysis. It uses modal coordinate equations
%    to compute the modal responses analytically, then convert the modal responses into
%    physical responses through the coordinate transformation.
%  Synopsis:
%    [eta,y,omega1,sdof2]=StepRespt(kk,mm,fd,u,t,C,q0,dq0,a,b)
%  Variable Description:
%    Input parameters:
%      kk, mm - System stiffness and mass matrices
%      fd - Input or forcing influence matrix
%      u - number of the loads
%      t - Time of evaluation
%      C - Output matrix
%      q0, dq0 - Initial conditions
%      a, b - Parameters for proportional damping [C]=a[M]+b[K]
%    Outpur parameters:
%      eta - modal coordinate response
%      y - physical coordinate response
%      omega - natrural frequency
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
V1=[V(:,1:sdof2)];                                    % truncate the modal vectors
%--------------------------------------------------------------------------
%  (3) Normalize the eigenvectors and compute parameters
%--------------------------------------------------------------------------
Factor=diag(V1'*mm*V1);
Vnorm=V1*inv(sqrt(diag(Factor)));                    %  eigenvectors are normalized
omega2=diag(sqrt(Vnorm'*kk*Vnorm));                        % natural frequencies
           
Fnorm=Vnorm'*fd;                                      % modal input force vector
%--------------------------------------------------------------------------
%  (4) Compute modal damping matrix from the proportional damping matrix
%--------------------------------------------------------------------------
Modamp=Vnorm'*(a*mm+b*kk)*Vnorm;                % form the Rayleigh damping
zeta=diag((1/2)*Modamp*inv(diag(omega2)));                     % the damping ratio

if (max(zeta) >= 1),
  disp('Warning - Your maximum damping ratio is grater than or equal to 1')
  disp('You have to reselect a and b ')
  pause
  disp('If you want to continue, type return key')
end
%--------------------------------------------------------------------------
%  (5) Find out step response of each modal coordinate analytically
%--------------------------------------------------------------------------
eta0=Vnorm'*mm*q0; deta0=Vnorm'*mm*dq0;
                % initial conditions for modal coordinates both displacement and velocity
eta=zeros(nstep,sdof2);

  for i=1:sdof2                               % responses are calculated for n modes
    omegad=omega(i)*sqrt(1-zeta(i)^2);
    phase=omegad*t;
    Exx=exp(-zeta(i)*omega(i)*t);

    C1=eta0(i);
    C2=(deta0(i)+eta0(i)*zeta(i)*omega(i))/omegad;

    D1=zeta(i)*omega(i)/omegad;
    II=ones(nstep,1);

    XX=Fnorm(i,u)/(omegad^2+zeta(i)^2*omega(i)^2);

    eta(:,i)=C1*Exx.*cos(phase)+C2*Exx.*sin(phase)...
             + XX*(II-Exx.*cos(phase)-D1*Exx.*sin(phase));
                                                % response for the step excitation
  end
%--------------------------------------------------------------------------
%  (6) Convert modal coordinate responses to physical coordinate responses
%--------------------------------------------------------------------------
eta=eta';
y=C*Vnorm*eta;

if (a+b)==0
  disp('The response results of undamping system')
else
  disp('The response results of damping system')
end

disp('The excitation is step force')
%--------------------------------------------------------------------------
%     The end
%--------------------------------------------------------------------------
