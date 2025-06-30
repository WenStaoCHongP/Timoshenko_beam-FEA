%-----------------------------------------------------------------------------------
% Ex.8.2.3
% to execute pole assignment design of frame structure
%-----------------------------------------------------------------------------------

M=4e5*eye(3); k1=2e8; k2=k1; k3=k2;
K=[k1+k2, -k2, 0; -k2, k2+k3, -k3; 0, -k3, k3];
C=0.7334*M+0.0026*K;
Bs=[1, -1, 0; 0, 1, -1; 0, 0, 1];
Ds=-M*[1;1;1];
A=[zeros(3), eye(3); -inv(M)*K, -inv(M)*C];
B=[zeros(3); inv(M)*Bs];
C0=eye(6);
D0=0;
D=[zeros(3,1);inv(M)*Ds];
lambda1=[-1.8218+9.9378*i,  -1.8218-9.9378*i,  -5.1045+27.8451*i,... 
-5.1045-27.8451*i   -7.5273+40.2094*i,  -7.5273-40.2094*i];
Ph1=inv(lambda1(1)*eye(6)-A)*B;
Ph2=inv(lambda1(2)*eye(6)-A)*B;
Ph3=inv(lambda1(3)*eye(6)-A)*B;  
Ph4=inv(lambda1(4)*eye(6)-A)*B;
Ph5=inv(lambda1(5)*eye(6)-A)*B;
Ph6=inv(lambda1(6)*eye(6)-A)*B;  
e=[1, 1, 0, 0, 0, 0; 0, 0, 1, 1, 0, 0; 0, 0, 0, 0, 1, 1];
Ke=-e*inv([Ph1(:,1)  Ph2(:,1)  Ph3(:,2)  Ph4(:,2)  Ph5(:,3)  Ph6(:,3)]);
Ke=real(Ke);
load xg.dat
xg=xg/(max(xg))*2;
X0=zeros(6,1);
G=ss(A,D, C0, D0);
G1=ss(A-B*Ke, D, C0, D0);
t=0:0.01:15;
[y, t, x]=lsim(G, xg, t, X0);
[y1, t, x1]=lsim(G1, xg, t, X0);
plot(t, x(:,1), '-.', t, x1(:,1))
legend('before control', 'after control')
xlabel('时间/s')
ylabel('第一层位移/cm',)
%-----------------------------------------------------------------------------------
