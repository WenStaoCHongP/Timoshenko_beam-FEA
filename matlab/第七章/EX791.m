E=2.1e11;      	%elastic molulus
poisson =0.3;         % poisson ratio
density=7.3e3;        %density
t=0.05;               %plate thickness
lx=2;                 %length in x direction
ly=2;                 %length in y direction
jdx=11;               %number of nodes in x direction
jdy=11;               %number of nodes in y direction

k(1:330,1:330)=0;     %system stiffness matrix
m(1:330,1:330)=0;    %system mass matrix

%prepare the arrays which are needed to describe this problem
en(1:100,1:4)=0;       %element node  
for ni=1:jdx-1
    for nj=1:jdy-1
        en(ni+(nj-1)*(jdx-1),1)=ni+(nj-1)*jdx;
        en(ni+(nj-1)*(jdx-1),2)=ni+1+(nj-1)*jdx;
        en(ni+(nj-1)*(jdx-1),4)=ni+nj*jdx;
        en(ni+(nj-1)*(jdx-1),3)=ni+1+nj*jdx;
    end
end
disp(1:jdx*jdy,1:3)=1;     % node displacement
constraints=1:jdx:jdx*jdy;  % constraints
disp(constraints,:)=0;
dof=0;                   %degree of freedom
for ni=1:jdx*jdy
    for nj=1:3
        if disp(ni,nj)~=0
            dof=dof+1;
            disp(ni,nj)=dof;
        end
    end
end

el=lx/(jdx-1);        %element length 
eh=ly/(jdy-1);       %element height
[ek,dm]=km(el/2,eh/2,mu,poisson,E,density); 
%km: function used to compute element stifness and mass matrix
%in this case, all elements have the same element stifness and mass matrix. 

%built system stifness and mass matrix.  
index(1:12)=0; % vector sontaining system dofs of nodes in each element. 
for loopi=1:(jdx-1)*(jdy-1)
    for zi=1:4
        index((zi-1)*3+1)=disp(en(loopi,zi),1);
        index((zi-1)*3+2)=disp(en(loopi,zi),2);
        index((zi-1)*3+3)=disp(en(loopi,zi),3);
    end
    for jx=1:12
        for jy=1:12
            if(index(jx)*index(jy)~=0)
                  k(index(jx),index(jy))=k(index(jx),index(jy))+ek(jx,jy);
                  m(index(jx),index(jy))=m(index(jx),index(jy))+em(jx,jy);
            end
        end
    end
end

%solve eigenvalue problem
[v,d] = eig(k,m);
tempd=diag(d);
[nd,sortindex]=sort(tempd);
v=v(:,sortindex);
mode_number=1:15;
frequency(mode_number)=sqrt(nd(mode_number))/(2*pi);

function [k,m]=km(a,b,poisson,t,E,density)
k=[E*t^3/(360-360*poisson^2)/a/b*(21-6*poisson+30*b^2/a^2+30*a^2/b^2),           E*t^3/(360-360*poisson^2)/a/b*(3*b+12*poisson*b+30*a^2/b),          E*t^3/(360-360*poisson^2)/a/b*(-3*a-12*poisson*a-30*b^2/a), E*t^3/(360-360*poisson^2)/a/b*(-21+6*poisson-30*b^2/a^2+15*a^2/b^2),          E*t^3/(360-360*poisson^2)/a/b*(-3*b-12*poisson*b+15*a^2/b),           E*t^3/(360-360*poisson^2)/a/b*(-3*a+3*poisson*a-30*b^2/a),  E*t^3/(360-360*poisson^2)/a/b*(21-6*poisson-15*b^2/a^2-15*a^2/b^2),           E*t^3/(360-360*poisson^2)/a/b*(-3*b+3*poisson*b+15*a^2/b),            E*t^3/(360-360*poisson^2)/a/b*(3*a-3*poisson*a-15*b^2/a), E*t^3/(360-360*poisson^2)/a/b*(-21+6*poisson+15*b^2/a^2-30*a^2/b^2),            E*t^3/(360-360*poisson^2)/a/b*(3*b-3*poisson*b+30*a^2/b),           E*t^3/(360-360*poisson^2)/a/b*(3*a+12*poisson*a-15*b^2/a);
E*t^3/(360-360*poisson^2)/a/b*(3*b+12*poisson*b+30*a^2/b),          E*t^3/(360-360*poisson^2)/a/b*(8*b^2-8*poisson*b^2+40*a^2),                               -30*E*t^3/(360-360*poisson^2)*poisson,          E*t^3/(360-360*poisson^2)/a/b*(-3*b-12*poisson*b+15*a^2/b),         E*t^3/(360-360*poisson^2)/a/b*(-8*b^2+8*poisson*b^2+20*a^2),                                                                 0,            
E*t^3/(360-360*poisson^2)/a/b*(3*b-3*poisson*b-15*a^2/b),          E*t^3/(360-360*poisson^2)/a/b*(2*b^2-2*poisson*b^2+10*a^2),                                                                 0,           
E*t^3/(360-360*poisson^2)/a/b*(-3*b+3*poisson*b-30*a^2/b),         E*t^3/(360-360*poisson^2)/a/b*(-2*b^2+2*poisson*b^2+20*a^2),                                                                 0;
E*t^3/(360-360*poisson^2)/a/b*(-3*a-12*poisson*a-30*b^2/a),                               -30*E*t^3/(360-360*poisson^2)*poisson,          E*t^3/(360-360*poisson^2)/a/b*(8*a^2-8*poisson*a^2+40*b^2),            E*t^3/(360-360*poisson^2)/a/b*(3*a-3*poisson*a+30*b^2/a),                                                                 0,         
E*t^3/(360-360*poisson^2)/a/b*(-2*a^2+2*poisson*a^2+20*b^2),           E*t^3/(360-360*poisson^2)/a/b*(-3*a+3*poisson*a+15*b^2/a),                                                                 0,          
E*t^3/(360-360*poisson^2)/a/b*(2*a^2-2*poisson*a^2+10*b^2),           E*t^3/(360-360*poisson^2)/a/b*(3*a+12*poisson*a-15*b^2/a),                                                                 0,         
E*t^3/(360-360*poisson^2)/a/b*(-8*a^2+8*poisson*a^2+20*b^2);
E*t^3/(360-360*poisson^2)/a/b*(-21+6*poisson-30*b^2/a^2+15*a^2/b^2),          E*t^3/(360-360*poisson^2)/a/b*(-3*b-12*poisson*b+15*a^2/b),            E*t^3/(360-360*poisson^2)/a/b*(3*a-3*poisson*a+30*b^2/a),  E*t^3/(360-360*poisson^2)/a/b*(21-6*poisson+30*b^2/a^2+30*a^2/b^2),           E*t^3/(360-360*poisson^2)/a/b*(3*b+12*poisson*b+30*a^2/b),           E*t^3/(360-360*poisson^2)/a/b*(3*a+12*poisson*a+30*b^2/a), E*t^3/(360-360*poisson^2)/a/b*(-21+6*poisson+15*b^2/a^2-30*a^2/b^2),            E*t^3/(360-360*poisson^2)/a/b*(3*b-3*poisson*b+30*a^2/b),          E*t^3/(360-360*poisson^2)/a/b*(-3*a-12*poisson*a+15*b^2/a),  E*t^3/(360-360*poisson^2)/a/b*(21-6*poisson-15*b^2/a^2-15*a^2/b^2),           E*t^3/(360-360*poisson^2)/a/b*(-3*b+3*poisson*b+15*a^2/b),           E*t^3/(360-360*poisson^2)/a/b*(-3*a+3*poisson*a+15*b^2/a);
E*t^3/(360-360*poisson^2)/a/b*(-3*b-12*poisson*b+15*a^2/b),         E*t^3/(360-360*poisson^2)/a/b*(-8*b^2+8*poisson*b^2+20*a^2),                                                                 0,
E*t^3/(360-360*poisson^2)/a/b*(3*b+12*poisson*b+30*a^2/b),          E*t^3/(360-360*poisson^2)/a/b*(8*b^2-8*poisson*b^2+40*a^2),                                30*E*t^3/(360-360*poisson^2)*poisson,           E*t^3/(360-360*poisson^2)/a/b*(-3*b+3*poisson*b-30*a^2/b),         E*t^3/(360-360*poisson^2)/a/b*(-2*b^2+2*poisson*b^2+20*a^2),                                                                 0,            
E*t^3/(360-360*poisson^2)/a/b*(3*b-3*poisson*b-15*a^2/b),          E*t^3/(360-360*poisson^2)/a/b*(2*b^2-2*poisson*b^2+10*a^2),                                                                 0;
E*t^3/(360-360*poisson^2)/a/b*(-3*a+3*poisson*a-30*b^2/a),                                                                 0,         
E*t^3/(360-360*poisson^2)/a/b*(-2*a^2+2*poisson*a^2+20*b^2),           E*t^3/(360-360*poisson^2)/a/b*(3*a+12*poisson*a+30*b^2/a),                                30*E*t^3/(360-360*poisson^2)*poisson,          E*t^3/(360-360*poisson^2)/a/b*(8*a^2-8*poisson*a^2+40*b^2),          E*t^3/(360-360*poisson^2)/a/b*(-3*a-12*poisson*a+15*b^2/a),                                                                 0,         
E*t^3/(360-360*poisson^2)/a/b*(-8*a^2+8*poisson*a^2+20*b^2),            E*t^3/(360-360*poisson^2)/a/b*(3*a-3*poisson*a-15*b^2/a),                                                                 0,          
E*t^3/(360-360*poisson^2)/a/b*(2*a^2-2*poisson*a^2+10*b^2);
E*t^3/(360-360*poisson^2)/a/b*(21-6*poisson-15*b^2/a^2-15*a^2/b^2),            E*t^3/(360-360*poisson^2)/a/b*(3*b-3*poisson*b-15*a^2/b),           E*t^3/(360-360*poisson^2)/a/b*(-3*a+3*poisson*a+15*b^2/a), E*t^3/(360-360*poisson^2)/a/b*(-21+6*poisson+15*b^2/a^2-30*a^2/b^2),           E*t^3/(360-360*poisson^2)/a/b*(-3*b+3*poisson*b-30*a^2/b),          E*t^3/(360-360*poisson^2)/a/b*(-3*a-12*poisson*a+15*b^2/a),  E*t^3/(360-360*poisson^2)/a/b*(21-6*poisson+30*b^2/a^2+30*a^2/b^2),          E*t^3/(360-360*poisson^2)/a/b*(-3*b-12*poisson*b-30*a^2/b),           E*t^3/(360-360*poisson^2)/a/b*(3*a+12*poisson*a+30*b^2/a), E*t^3/(360-360*poisson^2)/a/b*(-21+6*poisson-30*b^2/a^2+15*a^2/b^2),           E*t^3/(360-360*poisson^2)/a/b*(3*b+12*poisson*b-15*a^2/b),            E*t^3/(360-360*poisson^2)/a/b*(3*a-3*poisson*a+30*b^2/a);
E*t^3/(360-360*poisson^2)/a/b*(-3*b+3*poisson*b+15*a^2/b),          E*t^3/(360-360*poisson^2)/a/b*(2*b^2-2*poisson*b^2+10*a^2),                                                                 0,            
E*t^3/(360-360*poisson^2)/a/b*(3*b-3*poisson*b+30*a^2/b),         E*t^3/(360-360*poisson^2)/a/b*(-2*b^2+2*poisson*b^2+20*a^2),                                                                 0,          
E*t^3/(360-360*poisson^2)/a/b*(-3*b-12*poisson*b-30*a^2/b),          E*t^3/(360-360*poisson^2)/a/b*(8*b^2-8*poisson*b^2+40*a^2),                               -30*E*t^3/(360-360*poisson^2)*poisson,           E*t^3/(360-360*poisson^2)/a/b*(3*b+12*poisson*b-15*a^2/b),         E*t^3/(360-360*poisson^2)/a/b*(-8*b^2+8*poisson*b^2+20*a^2),                                                                 0;
E*t^3/(360-360*poisson^2)/a/b*(3*a-3*poisson*a-15*b^2/a),                                                                 0,          
E*t^3/(360-360*poisson^2)/a/b*(2*a^2-2*poisson*a^2+10*b^2),          E*t^3/(360-360*poisson^2)/a/b*(-3*a-12*poisson*a+15*b^2/a),                                                                 0,         
E*t^3/(360-360*poisson^2)/a/b*(-8*a^2+8*poisson*a^2+20*b^2),           E*t^3/(360-360*poisson^2)/a/b*(3*a+12*poisson*a+30*b^2/a),                               -30*E*t^3/(360-360*poisson^2)*poisson,          E*t^3/(360-360*poisson^2)/a/b*(8*a^2-8*poisson*a^2+40*b^2),           E*t^3/(360-360*poisson^2)/a/b*(-3*a+3*poisson*a-30*b^2/a),                                                                 0,         
E*t^3/(360-360*poisson^2)/a/b*(-2*a^2+2*poisson*a^2+20*b^2);
E*t^3/(360-360*poisson^2)/a/b*(-21+6*poisson+15*b^2/a^2-30*a^2/b^2),           E*t^3/(360-360*poisson^2)/a/b*(-3*b+3*poisson*b-30*a^2/b),           E*t^3/(360-360*poisson^2)/a/b*(3*a+12*poisson*a-15*b^2/a),  E*t^3/(360-360*poisson^2)/a/b*(21-6*poisson-15*b^2/a^2-15*a^2/b^2),            E*t^3/(360-360*poisson^2)/a/b*(3*b-3*poisson*b-15*a^2/b),            E*t^3/(360-360*poisson^2)/a/b*(3*a-3*poisson*a-15*b^2/a), E*t^3/(360-360*poisson^2)/a/b*(-21+6*poisson-30*b^2/a^2+15*a^2/b^2),           E*t^3/(360-360*poisson^2)/a/b*(3*b+12*poisson*b-15*a^2/b),           E*t^3/(360-360*poisson^2)/a/b*(-3*a+3*poisson*a-30*b^2/a),  E*t^3/(360-360*poisson^2)/a/b*(21-6*poisson+30*b^2/a^2+30*a^2/b^2),          E*t^3/(360-360*poisson^2)/a/b*(-3*b-12*poisson*b-30*a^2/b),          E*t^3/(360-360*poisson^2)/a/b*(-3*a-12*poisson*a-30*b^2/a);
E*t^3/(360-360*poisson^2)/a/b*(3*b-3*poisson*b+30*a^2/b),         E*t^3/(360-360*poisson^2)/a/b*(-2*b^2+2*poisson*b^2+20*a^2),                                                                 0,           
E*t^3/(360-360*poisson^2)/a/b*(-3*b+3*poisson*b+15*a^2/b),          E*t^3/(360-360*poisson^2)/a/b*(2*b^2-2*poisson*b^2+10*a^2),                                                                 0,           
E*t^3/(360-360*poisson^2)/a/b*(3*b+12*poisson*b-15*a^2/b),         E*t^3/(360-360*poisson^2)/a/b*(-8*b^2+8*poisson*b^2+20*a^2),                                                                 0,          
E*t^3/(360-360*poisson^2)/a/b*(-3*b-12*poisson*b-30*a^2/b),          E*t^3/(360-360*poisson^2)/a/b*(8*b^2-8*poisson*b^2+40*a^2),                                30*E*t^3/(360-360*poisson^2)*poisson;
E*t^3/(360-360*poisson^2)/a/b*(3*a+12*poisson*a-15*b^2/a),                                                                 0,         
E*t^3/(360-360*poisson^2)/a/b*(-8*a^2+8*poisson*a^2+20*b^2),           E*t^3/(360-360*poisson^2)/a/b*(-3*a+3*poisson*a+15*b^2/a),                                                                 0,          
E*t^3/(360-360*poisson^2)/a/b*(2*a^2-2*poisson*a^2+10*b^2),            E*t^3/(360-360*poisson^2)/a/b*(3*a-3*poisson*a+30*b^2/a),                                                                 0,         
E*t^3/(360-360*poisson^2)/a/b*(-2*a^2+2*poisson*a^2+20*b^2),          E*t^3/(360-360*poisson^2)/a/b*(-3*a-12*poisson*a-30*b^2/a),                                30*E*t^3/(360-360*poisson^2)*poisson,          E*t^3/(360-360*poisson^2)/a/b*(8*a^2-8*poisson*a^2+40*b^2)];
w=a*b*t*density;
syms kx yt kxi yti real;
ni=1/8*(1+kx*kxi)*(1+yt*yti)*(2+kx*kxi+yt*yti-kx^2-yt^2);
nix=-1/8*b*yti*(1+kx*kxi)*(1+yt*yti)*(1-yt^2);
niy=1/8*a*kxi*(1+kx*kxi)*(1+yt*yti)*(1-kx^2);
n(1)=subs(ni,{kxi,yti},{-1,-1});
n(2)=subs(nix,{kxi,yti},{-1,-1});
n(3)=subs(niy,{kxi,yti},{-1,-1});
 
n(4)=subs(ni,{kxi,yti},{1,-1});
n(5)=subs(nix,{kxi,yti},{1,-1});
n(6)=subs(niy,{kxi,yti},{1,-1});
 
n(7)=subs(ni,{kxi,yti},{1,1});
n(8)=subs(nix,{kxi,yti},{1,1});
n(9)=subs(niy,{kxi,yti},{1,1});
 
n(10)=subs(ni,{kxi,yti},{-1,1});
n(11)=subs(nix,{kxi,yti},{-1,1});
n(12)=subs(niy,{kxi,yti},{-1,1});
 
temp=n'*n;
m1=int(temp,kx,-1,1);
m=int(m1,yt,-1,1);
m=m*w;
m=double(m);
