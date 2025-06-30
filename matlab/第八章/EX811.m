%-------------------------------------------------------------------------------------% Example 8.1.1 
% to seek the maximal value of the first natural frequency for plane 2-bar truss
%-------------------------------------------------------------------------------------
 A0=[5 8];                               % initial cross-section area
v=ones(2,1);
k=1;
flag=1;
D(1)=0;
B0=5;
while flag
   m=3.12*A0(1)+3.9*A0(2);
M=[m 0;0 m];
k1=375*A0(1)+192*A0(2);
k2=144*A0(2);
k3=108*A0(2);
K=7.0e5*[k1 k2;k2 k3];
[V,D1]=eig(K,M);                                % compute eigenvalue and eigenvector 
[lambda,k4]=sort(diag(real((D1))));
V=real(V(:,k4));
Factor=diag(V'*M*V);
Vnorm=V*inv(sqrt(diag(Factor)));                   % normalize eigenvector
freq=sqrt(lambda)/(2*pi);                  
A0=[1  1];
for i=1:2
m=3.12*A0(1)+3.9*A0(2);
M=[m 0;0 m];
k1=375*A0(1)+192*A0(2);
k2=144*A0(2);
k3=108*A0(2);
K=7.0e5*[k1 k2;k2 k3];
v(i)=Vnorm(:,1)'*(K-lambda(1)*M)*Vnorm(:,1);     % sensitivity analysis
end
   sum=0.0;
   for i=1:2
      sum=sum+v(i)^2;
   end
   sum=sqrt(sum);
   for i=1:2
      v(i)=-v(i)/sum;
      n=0;
      if v(i)<0
         n=n+1;
         c0(n)=A0(i)/v(i);
     end
   end
     
c0max=max(c0);
B=c0max/B0; 
 
for i=1:2
   A(i)=A0(i)+B*v(i);
end
A0=A;
v0=v;
k=k+1;
D(k)=lambda(1)
epilson=(abs(D(k)-D(k-1)))/D(k);

if epilson<=1e-5
   flag=0;
else
   flag=1;
end
end

plot(D(2:k),'--*')
xlabel(' 迭代次数')
ylabel('第一特征值(/s^2)')
