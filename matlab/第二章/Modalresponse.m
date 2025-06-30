function  varargout=Modalresponse(varargin)
%-------------------------------------------------------------------------
% Purpose:
% This function is used to calculate the modal and physical response for a
% strcutral system.  
% Synopsis:
%       [eta,yim]=Modalresponse(M,K,Fu,u,t,Cy,q0,dq0,nummode)
%       [eta,yim]=Modalresponse(M,K,Fu,u,t,Cy,q0,dq0,a,b,nummode)       
% Variable Description:
%       Input parameters
%       K - System siffness matrix
%       M - Sytem mass matrix
%       Fu - Force influence matrix
%       u  - Force vector
%       t  - Time of evaluation
%       Cy - Output matrix
%       q0,dq0 - Initial conditions
%       nummode - the number of extracted modes
%       a, b - Parameters for proportional damping [C]=a[M]+b[K]
%       Output parameters
%       eta - modal coordinate and velocity response
%       yim  -  physcial coordinate response
%    
% Author Dr.XU Bin, Time:2006-11-29
%--------------------------------------------------------------------------
disp('')
disp('Please wait!! - The job is being performed.')
%-------------------------------------------------------------------------
 if nargin<9 |nargin>11|nargin==10
    error('Incorrect number of input arguments')
 else
     switch nargin
         case 9
 %------------------------------------------------------------------------
 % Solve the eigenvalue problem and normalized the eigenvectors
 %-----------------------------------------------------------------------
  M=varargin{1};
  K=varargin{2};
  Fu=varargin{3};
  t=varargin{5};
  [n,n]=size(M);[n,m]=size(Fu);
  nstep=size(t');
  [V,D]=eig(K,M);
  [lambda,k]=sort(diag(D));      % Sort the eigenvalues and eigenvalues                      
  V=V(:,k);
  Factor=diag(V'*M*V);
  Vnorm=V*inv(sqrt(diag(Factor))); % Eigenvectors are normalized
  q0=varargin{7};               
  dq0=varargin{8};               
  nummode=varargin{9};
  eta0=Vnorm'*q0;              % Initial conditions for modal coordinates
  deta0=Vnorm'*dq0;            % both displacement and velocity
  eta=zeros(nstep,nummode);
  for i=1:nummode
          u=varargin{4};
          A=[0 1;-lambda(i) 0];
          B=[zeros(1,m);Vnorm(:,i)'*Fu];
          C=eye(2);
          D=0;
          x0=[eta0(i);deta0(i)];
          x=lsim(ss(A,B,C,D),u,t,x0);
          for j=1:nstep
              eta(j,i)=x(j,1);
          end
  end
  Cy=varargin{6};
  yim=Cy*Vnorm(:,1:nummode)*eta'; % Convert modal coordinate response to physical coodrinate responses
  varargout{1}=eta;
  varargout{2}=yim;
         
         case 11
 %------------------------------------------------------------------------
 % Solve the eigenvalue problem and normalized the eigenvectors
 %-----------------------------------------------------------------------
  M=varargin{1};
  k=varargin{2};
  Fu=varargin{3};
  t=varargin{5};
  [n,n]=size(M);[n,m]=size(Fu);
  nstep=size(t');
  [V,D]=eig(K,M);
  [lambda,k]=sort(diag(D));      % Sort the eigenvalues and eigenvalues                      
  V=V(:,k);
  Factor=diag(V'*M*V);
  Vnorm=V*inv(sqrt(diag(Factor))); % Eigenvectors are normalized
  Omega=diag(sqrt(Vnorm'*K*Vnorm));% Natural frequencies
 %------------------------------------------------------------------
 % Compute modal damping matrix from the proportional damping matrix
 %------------------------------------------------------------------
  a=varargin{9};
  b=varargin{10};
  Modamp=Vnorm'*(a*M+b*K)*Vnorm;
  zeta=diag((1/2)*Modamp*inv(diag(Omega)));
  if (max(zeta)>=1)
      disp( 'Warning - Maximum damping ratio is greater than or equal to 1')
      disp('You have to reselect a and b')
      pause
      disp('If you want to continue, type return key')
  end
  %-----------------------------------------------------------------------
  q0=varargin{7};               
  dq0=varargin{8};               
  nummode=varargin{11};
  eta0=Vnorm'*M*q0;              % Initial conditions for modal coordinates
  deta0=Vnorm'*M*dq0;            % - both displacement and velocity
  eta=zeros(nstep,nummode);
  for i=1:nummode
  
          u=varargin{4};
          A=[0 1;-lambda(i) -Modamp(i)];
          B=[zeros(1,m);Vnorm(:,i)'*Fu];
          C=eye(2);
          D=0;
          x0=[eta0(i);deta0(i)];
          x=lsim(ss(A,B,C,D),u,t,x0);
          for j=1:nstep
              eta(j,i)=x(j,1);
          end
  end
  Cy=varargin{6};
  yim=Cy*Vnorm(:,1:nummode)*eta'; % Convert modal coordinate response to physical coodrinate %responses
  varargout{1}=eta;
  varargout{2}=yim;
     end
 end
 %---------------------------------------------------------------------------------------------------------------------------------------
