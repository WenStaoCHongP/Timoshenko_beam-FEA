%-----------------------------------------------------------------------------------
% Ex.8.2.2
% to execute LQG and LQR designs of frame structure
%-----------------------------------------------------------------------------------
M=4e5*eye(3); k1=2e8; k2=k1; k3=k2;
K=[k1+k2, -k2, 0; -k2, k2+k3, -k3; 0, -k3, k3];
C=0.7334*M+0.0026*K;
Bs=[1, -1, 0; 0, 1, -1; 0, 0, 1];
Ds=-M*[1;1;1];
A=[zeros(3), eye(3); -inv(M)*K, -inv(M)*C];
B=[zeros(3); inv(M)*Bs];
D=[zeros(3,1);inv(M)*Ds];
Q=100*[K, zeros(3); zeros(3), M];
R=5e-6*eye(3);
G=lqr(A, B, Q, R);
C0=[0 0 0 1 0 0;0 0 0 0 1 0;0 0 0 0 0 1];
Qe=1e-4;
Re=1e-2*eye(3);
Ke=lqe2(A, D, C0, Qe, Re);
Am=A-B*G-Ke*C0*A+Ke*C0*B*G;
Amm=inv(eye(6)-Ke*C0)*Am;
Bmm=D;
load Xg.dat
Xg=Xg/max(Xg)*2;
X0=zeros(6,1);
C1=eye(6);
D0=0;
G1=ss(Amm,Bmm, C1, D0);
t=0:0.02:30;
[y, t, x]=lsim(G1, Xg, t, X0);
G2=ss(A-B*G,Bmm,C1,D0);
[y1, t1, x1]=lsim(G2, Xg, t, X0);
hold on
plot(t, x(:,1), '*')
plot(t, x1(:,1))
hold off
legend('LQG control', 'LQR control')
xlabel('时间/s')
ylabel('第一层位移/m')
