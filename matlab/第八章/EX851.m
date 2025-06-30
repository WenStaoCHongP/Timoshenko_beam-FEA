E=2*10^11;						
L=1;								
m1=100;						
R=0.03;								
e=2*10^(-4);							
g1=9.8;								
w=150;								
dt=pi/500;								
t0=0;									
I=pi*R^4/4;							
kz=48*E*I/L^3;							
wn=(kz/m1)^0.5;							
q=m1*g1/kz;							
 
U=e/q;								
W=w/wn;								
 
r=0.5;
a=0.25;
a1=1/(a*dt^2);
a2=-1/(a*dt);
a3=-(0.5/a-1);
b1=r/(a*dt);
b2=1-r/a;
b3=(1-0.5*r/a)*dt;						
 
m=eye(2);								
c=2*[0.01/W 0;0 0.01/W];					
k=[1/W^2 0;0 1/W^2];						
g=[U*sin(t0);U*cos(t0)];					
 
u0=[0;0];
up0=[0.001;0.001];
upp0=m^(-1)*(g-k*u0-c*up0);				
 
u=u0;
up=up0;
upp=upp0;
fid=fopen('zhuanzi.txt','w');					
for t=(t0+dt):dt:300*pi						
    g=[U*sin(t);U*cos(t)];					
    K=a1*m+b1*c+k;
    G=g+m*(a1*u-a2*up-a3*upp)+c*(b1*u-b2*up-b3*upp);
    u1=K^(-1)*G;
    up1=b1*(u1-u)+b2*up+b3*upp;
    upp1=a1*(u1-u)+a2*up+a3*upp;			
    u=u1;
    up=up1;
    upp=upp1;
    fprintf(fid,'%f     %f      %f\n',u1(1),u1(2),t);		
end
fclose(fid);				
