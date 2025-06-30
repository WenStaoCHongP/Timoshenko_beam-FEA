function varargout=freqresponse(varargin)
%-------------------------------------------------------------------------
% Purpose:
% This function is used to calculate the freqency response for a strcutral 
% system.  
% Synopsis:
%       Hd=freqresponse(K,M,C,Omega,nummode,'displacement')
%       Hv=freqresponse(K,M,C,Omega,nummode,'velocity')
%       Ha=freqresponse(K,M,C,Omega,nummode,'acceleration')
%       [Hd,Hv,Ha]=freqresponse(K,M,C,Omega,nummode,'all')
%       Hd=freqresponse(K,M,C,Omega,nummode,'displacement',p,q)
%       Hv=freqresponse(K,M,C,Omega,nummode,'velocity',p,q)
%       Ha=freqresponse(K,M,C,Omega,nummode,'acceleration',p,q)
%       [Hd,Hv,Ha]=freqresponse(K,M,C,Omega,nummode,'all',p,q)
% Variable Description:
%       Input parameters
%       K - System siffness matrix
%       M - Sytem mass matrix
%       C - System damping matrix
%       Omega - Frequency range
%       nummode - the number of extracted modes
%       p/q    -- Serial number of the degrees of freedom
%       Output parameters
%       Hd - Displacement freqency response
%       Hv - Velocity freqency response
%       Ha - Acceleration freqency response
% Author Dr.XU Bin, Time:2006-11-29
%--------------------------------------------------------------------------

if nargin<6 |nargin>8|nargin==7
    error('Incorrect number of input arguments')
else
K=varargin{1};
M=varargin{2};
C=varargin{3};
Omega=varargin{4};
nummode=varargin{5};
nsdof=length(diag(K));
[V,D]=eig(K,M);
[lambda,k]=sort(diag(D))
V=V(:,k);
Factor=diag(V'*M*V);
Vnorm=V*inv(sqrt(diag(Factor)));
Dampr=Vnorm'*C*Vnorm;
n=length(Omega);
switch nargin
   
case 6
 switch varargin{6}
     case 'displacement'
    Hd1=zeros(nsdod,nsdof);    
    Hd=zeros(nsdof,nsdof,n);
    for p=1:n
        for q=1:nummode
            Hd1= Hd1+Vnorm(:,q)*Vnorm(:,q)'/(lambda(q)-Omega(p)^2+i*Dampr(q,q)*Omega(p));
        end
        for l=1:nsdof
            for h=1:nsdof
                Hd(l,h,p)=Hd(l,h,p)+Hd1(l,h);
            end
        end         
    end
    varargout{1}=Hd;
     case 'velocity'
      Hv1=zeros(nsdof,nsdof);
      Hv=zeros(nsdof,nsdof,n);
      for p=1:n
           for q=1:nummode
             Hv1= Hv1+i*Omega(p)*Vnorm(:,q)*Vnorm(:,q)'/(lambda(q)-Omega(p)^2+i*Dampr(q,q)*Omega(p));
           end
           for l=1:nsdof
              for h=1:nsdof
                Hv(l,h,p)=Hv(l,h,p)+Hv1(l,h);
              end
          end         
      end
      varargout{1}=Hv;
     case 'acceleration'
       Ha1=zeros(nsdof,nsdof);
       Ha=zeros(nsdof,nsdof,n);
       for p=1:n
           for q=1:nummode
               Ha1= Ha1-Omega(p)^2*Vnorm(:,q)*Vnorm(:,q)'/(lambda(q)-Omega(p)^2+i*Dampr(q,q)*Omega(p));
           end
           for l=1:nsdof
             for h=1:nsdof
                Ha(l,h,p)=Ha(l,h,p)+Ha1(l,h);
             end
           end         
       end
       varargout{1}=Ha;
     case 'all'
            Hd1=zeros(nsdof,nsdof);
            Hv1=zeros(nsdof,nsdof);
            Ha1=zeros(nsdof,nsdof);
            Hd=zeros(nsdof,nsdof,n);
            Hv=zeros(nsdof,nsdof,n);
            Ha=zeros(nsdof,nsdof,n);
        for p=1:n
           for q=1:nummode
              
               Hd1= Hd1+Vnorm(:,q)*Vnorm(:,q)'/(lambda(q)-Omega(p)^2+i*Dampr(q,q)*Omega(p)); 
               Hv1= Hv1+i*Omega(p)*Vnorm(:,q)*Vnorm(:,q)'/(lambda(q)-Omega(p)^2+i*Dampr(q,q)*Omega(p));
               Ha1= Ha1-Omega(p)^2*Vnorm(:,q)*Vnorm(:,q)'/(lambda(q)-Omega(p)^2+i*Dampr(q,q)*Omega(p));
           end
           for l=1:nsdof
             for h=1:nsdof
                Hd(l,h,p)=Hd(l,h,p)+Hd1(l,h);
                Hv(l,h,p)=Hv(l,h,p)+Hv1(l,h);
                Ha(l,h,p)=Ha(l,h,p)+Ha1(l,h);
             end
           end         
        end
        varargout{1}=Hd;
        varargout{2}=Hv;
        varargout{3}=Ha;
 end
 
 
    case 8
      nodei=varargin{7};
      nodej=varargin{8};
    switch varargin{6}
        
     case 'displacement'
           Hd=zeros(n,1);
       for p=1:n
            for q=1:nummode
                 Hd(p)= Hd(p)+Vnorm(nodei,q)*Vnorm(nodej,q)'/(lambda(q)-Omega(p)^2+i*Dampr(q,q)*Omega(p));
            end      
       end
           varargout{1}=Hd;
     case 'velocity'
          Hv=zeros(n,1);
        for p=1:n
           for q=1:nummode
                Hv(p)= Hv(p)+i*Omega(p)*Vnorm(nodei,q)*Vnorm(nodej,q)'/(lambda(q)-Omega(p)^2+i*Dampr(q,q)*Omega(p));
           end
        end
           varargout{1}=Hv;
           
     case 'acceleration'
          Ha=zeros(n,1);
          for p=1:n
              for q=1:nummode
                  Ha(p)= Ha(p)-Omega(p)^2*Vnorm(nodei,q)*Vnorm(nodej,q)'/(lambda(q)-Omega(p)^2+i*Dampr(q,q)*Omega(p));
              end          
          end
          varargout{1}=Ha;
     case 'all'
           
            Hd=zeros(n,1);
            Hv=zeros(n,1);
            Ha=zeros(n,1);
        for p=1:n
           for q=1:nummode
        
               Hd(p)= Hd(p)+Vnorm(nodei,q)*Vnorm(nodej,q)'/(lambda(q)-Omega(p)^2+i*Dampr(q,q)*Omega(p)); 
               Hv(p)= Hv(p)+i*Omega(p)*Vnorm(nodei,q)*Vnorm(nodej,q)'/(lambda(q)-Omega(p)^2+i*Dampr(q,q)*Omega(p));
               Ha(p)= Ha(p)-Omega(p)^2*Vnorm(nodei,q)*Vnorm(nodej,q)'/(lambda(q)-Omega(p)^2+i*Dampr(q,q)*Omega(p));
           end
                  
        end
        varargout{1}=Hd;
        varargout{2}=Hv;
        varargout{3}=Ha;  
        
    end
end
end
%-------------------------------------------------------------------------------------

