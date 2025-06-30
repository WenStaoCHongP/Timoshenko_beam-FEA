function [v, d]=Inviter(K, M, p, epsilon)
%--------------------------------------------------------------------------
%  Purpose:
%     The function calculates the first p eigenvalues and eigenvectors for 
%     a structural system using inverse interation method.
%  Synopsis:
%     [v, d]=Inviter(K, M, p, epsilon)
%  Variable Description:
%     Input parameters
%     K - System siffness matrix
%     M - Sytem mass matrix
%     p - the number of calculated eigenvalues/eigenvectors
%     epsilon - the required precise
%     Output parameters
%     V - First p eigenvecotrs
%     d - First p eigenvalues
%--------------------------------------------------------------------------
[n1,n2]=size(K);
KM=zeros(n1,1);
for i=1:n1
    KM(i)=K(i,i)/M(i,i);
end
[Y,I]=sort(KM,1,'ascend');
y=zeros(n1,p);
y(:,1)=diag(M);
for i=2:p
    y(I(i-1),i)=1;
end
d=zeros(p,1);
v=zeros(n1,p);

   for j=1:p
        flag=1;
    while flag
        if j>1
           for i=1:(j-1)
               xs=v(:,i)'*y(:,j)*v(:,i);
               xs=M*xs;
              y(:,j)=y(:,j)-xs;
           end
        end
          
        xiter=K\y(:,j);
        yiter=M*xiter;
        diter=transpose(xiter)*y(:,j)/(transpose(xiter)*yiter);
        if (abs(diter-d(j))/diter)<epsilon  
            flag=0;
        else
            flag=1;
        end
           y(:,j)=yiter/((transpose(xiter)*yiter)^(1/2));
           d(j)=diter;     
      end
         d(j)=diter;
         v(:,j)=xiter/((transpose(xiter)*yiter)^(1/2));  
   end      
%---------------------------------------------------------------
%   The end 
%---------------------------------------------------------------     
