function [U,U2,G,freq,lambda]=FE(nelx,nely,xold)
%-------------------------------------------------------------------
% FE-ANALYSIS
%-------------------------------------------------------------------
[KE,Me] = lk; 
K = sparse(3*(nelx+1)*(nely+1), 3*(nelx+1)*(nely+1)); M= sparse(3*(nelx+1)*(nely+1), 3*(nelx+1)*(nely+1));
F = sparse(3*(nely+1)*(nelx+1),1); F1=sparse(3*(nely+1)*(nelx+1),3*(nely+1)*(nelx+1));
U = sparse(3*(nely+1)*(nelx+1),1);  U1 = sparse(3*(nely+1)*(nelx+1),3*(nely+1)*(nelx+1)); 
for elx = 1:nelx
  for ely = 1:nely
    n1 = (nely+1)*(elx-1)+ely; 
    n2 = (nely+1)* elx   +ely;
    edof = [3*n1-2;3*n1-1;3*n1;3*n2-2;3*n2-1;3*n2;3*n2+1;3*n2+2;3*n2+3;3*n1+1;3*n1+2;3*n1+3];
    K(edof,edof) = K(edof,edof) + ((xold(ely,elx)^3)/((0.01)^3))*KE;
    M(edof,edof) =M(edof,edof)+((xold(ely,elx))/(0.01))*Me;
  end
end
% DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
F(3*(0.5*nelx*(nely+1)+0.5*nely)+1,1) = -10000;
fixeddofs   = union([1,3*(nely+1)-2],[3*(nelx)*(nely+1)+1,3*(nelx+1)*(nely+1)-2]);
alldofs     = [1:3*(nely+1)*(nelx+1)];
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
V1=zeros(3*(nely+1)*(nelx+1),1);
V1(freedofs,:)=Vnorm(:,1);
V1(fixeddofs,:)=0;
for elx = 1:nelx
  for ely = 1:nely
    n1 = (nely+1)*(elx-1)+ely; 
    n2 = (nely+1)* elx   +ely;
    edof = [3*n1-2;3*n1-1;3*n1;3*n2-2;3*n2-1;3*n2;3*n2+1;3*n2+2;3*n2+3;3*n1+1;3*n1+2;3*n1+3];
    K1=3*((xold(ely,elx)^2)/((0.01)^3))*KE;
    M1=(1/(0.01))*Me;
    G(ely,elx)=[V1(edof,1)'*(K1-lambda(1)*M1)*V1(edof,1)];
   end
end
%--------------------------------------------------------------------------
for i=1:3*(nely+1)*(nelx+1)
     if((i==3*(0.5*nelx*(nely+1)+0.5*nely)-2)|(i==3*(0.5*nelx*(nely+1)+0.5*nely)+1))
     F1(i,i)=-1;
     U2([freedofs],i)=K(freedofs,freedofs) \ F1([freedofs],i);                                         % obtain unit virtual load vectors 
     U2([fixeddofs],i)= 0;
    end
   %  U1(:,[fixeddofs])=0;
end
