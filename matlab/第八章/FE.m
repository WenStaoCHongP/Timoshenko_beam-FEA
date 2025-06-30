function [U,U2,G,freq,lambda]=FE(nelx,nely,xold)
%------------------------------------------------------------
% FE-ANALYSIS
%------------------------------------------------------------
[KE,ME] = lk; 
K = sparse(2*(nelx+1)*(nely+1), 2*(nelx+1)*(nely+1)); M = sparse(2*(nelx+1)*(nely+1), 2*(nelx+1)*(nely+1)); 
F = sparse(2*(nely+1)*(nelx+1),1); F1=sparse(2*(nely+1)*(nelx+1),2*(nely+1)*(nelx+1));
U = sparse(2*(nely+1)*(nelx+1),1);  U2 = sparse(2*(nely+1)*(nelx+1),2*(nely+1)*(nelx+1)); 
for elx = 1:nelx
  for ely = 1:nely
    n1 = (nely+1)*(elx-1)+ely; 
    n2 = (nely+1)* elx   +ely;
    edof = [2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n2+1; 2*n2+2; 2*n1+1; 2*n1+2];
    K(edof,edof) = K(edof,edof) + xold(ely,elx)*KE;
    M(edof,edof) = M(edof,edof) + xold(ely,elx)*ME;
  end
end

% DEFINE LOADS AND SUPPORTS 
F(2*(nelx+0.5)*(nely+1)+1,1) = -5200;
F(2*(nelx+0.5)*(nely+1)-1,1)=-5200;
F(2*(nelx+0.5)*(nely+1)+3,1)=-5200;
fixeddofs   = [1:2*(nely+1)];
alldofs     = [1:2*(nely+1)*(nelx+1)];
freedofs    = setdiff(alldofs,fixeddofs);
% SOLVING
U(freedofs,:) = K(freedofs,freedofs) \ F(freedofs,:);  
U(fixeddofs,:)= 0;
% compute eigenvalues and eigenvectors
[V,D]=eigs(K(freedofs,freedofs), M(freedofs,freedofs),1,'SM');
[lambda,Kt]=sort(diag(D));
Factor=diag(V(:,1)'*M(freedofs,freedofs)*V(:,1));
Vnorm=V(:,1)*inv(sqrt(diag(Factor)));                         % normalize eigenvector
freq=sqrt(lambda)/(2*pi);
V1=zeros(2*(nely+1)*(nelx+1),1);
V1(freedofs,:)=Vnorm(:,1);                                   
V1(fixeddofs,:)=0;                                           % delete of zero element
% computation of frequency sensitivity
for elx = 1:nelx
  for ely = 1:nely
    n1 = (nely+1)*(elx-1)+ely; 
    n2 = (nely+1)* elx   +ely;
    edof = [2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n2+1; 2*n2+2; 2*n1+1; 2*n1+2];
    K1=KE;
    M1=ME;
    G(ely,elx)=[V1(edof,1)'*(K1-lambda(1)*M1)*V1(edof,1)];
   end
end
% compute unit virtual load vector
%-----------------------------------------------------------------
for i=1:2*(nely+1)*(nelx+1)
     if(i==2*(nelx+0.5)*(nely+1)+1)
     F1(i,i)=-1;
     U2([freedofs],i)=K(freedofs,freedofs) \ F1([freedofs],i);
     U2([fixeddofs],i)= 0;
     end
end
