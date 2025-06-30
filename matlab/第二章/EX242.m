%--------------------------------------------------------
%  Ex.2.2 plane 9-bar truss
%--------------------------------------------------------
%  To compute the dynamic response of node 4 on the vertical direction
%
%   
E=2.0e11;
A=2.5e-3;
density=7860; 
node_number=6;
elment_number=9;
nc=[0 0;4 0;4 3;8 0;8 3;12 0];            %node_coordinate
en=[1 2;1 3;2 3;2 4;3 4;3 5;4 5;4 6;5 6];   % element_node  
ed(1:node_number,1:2)=1;            % elment_displacement 
dof=0;
for loopi=1:node_number
     for loopj=1:2
            dof=dof+1;
            ed(loopi,loopj)=dof;
    end
end
ek=E*A*[1 0 -1 0;0 0 0 0;-1 0 1 0;0 0 0 0]; % elment_stiffness_matrix  
em=(density*A)/6*[2 0 1 0; ...           %elment_mass_matrix  
                   0 2 0 1; ...
                  1 0 2 0; ...
                  0 1 0 2];   
k(1:dof,1:dof)=0;                     %tructural_stiffness_matrix  
m=k;                                % structural_mass_matrix, same size  with  k 
theta(1:9)=0;
el(1:9)=0;
e2s(1:4)=0; %index of transform the elment displament number to structural
for loopi=1:elment_number
    for zi=1:2
        e2s((zi-1)*2+1)=ed(en(loopi,zi),1);
        e2s((zi-1)*2+2)=ed(en(loopi,zi),2);
    end
    el(loopi)=sqrt((nc(en(loopi,1),1)-nc(en(loopi,2),1))^2+(nc(en(loopi,1),2)-nc(en(loopi,2),2))^2);
    theta(loopi)=asin((nc(en(loopi,1),2)-nc(en(loopi,2),2))/el(loopi));
    lmd=[cos(theta(loopi)) sin(theta(loopi)); -sin(theta(loopi)) cos(theta(loopi))]; 
    t=[lmd zeros(2); zeros(2) lmd];
    dk=t'*ek*t/el(loopi);
    dm=t'*em*t*el(loopi);
    for jx=1:4
        for jy=1:4
                k(e2s(jx),e2s(jy))=k(e2s(jx),e2s(jy))+dk(jx,jy);
                m(e2s(jx),e2s(jy))=m(e2s(jx),e2s(jy))+dm(jx,jy);
        end
    end
end
c=0*k+0*m;
nt=1500;dt=0.0001;
time=0:dt:nt*dt;

q0=zeros(dof,1);
dq0=zeros(dof,1);
bcdof=zeros(dof,1);
fd=zeros(dof,nt);
for i=1:nt
    fd(6,i)=200;
end
[acc,vel,dsp]=TransResp1(k,c,m,fd,bcdof,nt,dt,q0,dq0);
%[acc,vel,dsp]=TransResp3(k,c,m,fd,bcdof,nt,dt,q0,dq0);
%[acc,vel,dsp]=TransResp4(k,c,m,fd,bcdof,nt,dt,q0,dq0);
%[acc,vel,dsp]=TransResp5(k,c,m,fd,bcdof,nt,dt,q0,dq0);
plot(time, dsp(8,:))
 xlabel('Time(seconds)')
 ylabel(' Vertical displ. (m)')

% state-space method
 AG=[zeros(dof) eye(dof);-inv(m)*k zeros(dof)];
 BG=[zeros(dof,1);
     inv(m)*[0 0 0 0 0 1 0 0 0]'];
 CG=eye(2*dof);
 DG=0;
nt=1500;
det=0.0001;
time=0:det:(nt*det);
for i=1:(nt+1)
      u(i)=200;
end
G=ss(AG,BG,CG,DG);
dsp=lsim(G,u,time);
plot(time,dsp(:,6))
xlabel('Time(seconds)')
ylabel('Tip displ. (m)')

% modal analysis method
nt=1500;
det=0.0001;
time=0:det:(nt*det);
for i=1:(nt+1)
      u(i)=200;
end
Fu=[0 0 0 0 0 1 0 0 0]';
Cy=eye(dof);
q0=zeros(dof,1);
dq0=zeros(dof,1);
nummode=9; 
[eta,dsp]=Modalresponse(m,k,Fu,u,time,Cy,q0,dq0,nummode);
plot(time,dsp(6,:))
xlabel('Time(seconds)')
ylabel(' Vertical displ. (m)')
