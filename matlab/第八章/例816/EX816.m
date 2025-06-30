%-------------------------------------------------------------------------------------
% Ex. 8.1.6
% four  angular point with hinged supported condition for square plate with multidisplacement constraints
%-------------------------------------------------------------------------------------
% INITIALIZE
nelx=30;
nely=18;
x(1:nely,1:nelx) = 0.05*ones(nely,nelx); 
a=0.1;
b=0.1;
p=7800;
loop = 0; 
h=0.5;
U4=-0.01;
change =0.02;
q=0.4;
U6=-0.0028;
d=1.4;
g=0.31;
freq1=3.5;
% START ITERATION
while change>0.0001
 loop = loop + 1;
  xold = x; 
% FE-ANALYSIS
  [U,U2,G,freq,lambda]=FE(nelx,nely,xold);
  if ((loop)>28),break,end 
% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
   [KE,Me] = lk; 
    c = 0.;
    U(3*(0.5*nelx*(nely+1)+0.5*nely)+1,1)
    U(3*(0.5*nelx*(nely+1)+0.5*nely)-2,1)
    freq(1)
  for ely = 1:nely
    for elx = 1:nelx
      n1 = (nely+1)*(elx-1)+ely; 
      n2 = (nely+1)* elx   +ely;
      edof = [3*n1-2;3*n1-1;3*n1;3*n2-2;3*n2-1;3*n2;3*n2+1;3*n2+2;3*n2+3;3*n1+1;3*n1+2;3*n1+3];
      c = c + 4*((xold(ely,elx)))*p*a*b;
      dc(ely,elx) = 4*p*a*b;
      Ue = U([3*n1-2;3*n1-1;3*n1;3*n2-2;3*n2-1;3*n2;3*n2+1;3*n2+2;3*n2+3;3*n1+1;3*n1+2;3*n1+3],1);
      F2(edof,:)=((xold(ely,elx)^3)/((0.01)^3))*KE*Ue;
      dU(ely,elx)=F2(edof,:)'*U2([edof],3*(0.5*nelx*(nely+1)+0.5*nely)+1);
      dU1(ely,elx)=F2(edof,:)'*U2([edof],3*(0.5*nelx*(nely+1)+0.5*nely)-2);
    end  
 end
   lamda1(loop)=c/(3*abs(U4));
   lamda2(loop)=c/(3*abs(U6));
   lamda1(loop+1)=lamda1(loop)*(U(3*(0.5*nelx*(nely+1)+0.5*nely)+1,1)/U4)^h;
   lamda2(loop+1)=lamda2(loop)*(U(3*(0.5*nelx*(nely+1)+0.5*nely)-2,1)/U6)^h;
   lamda5=lamda1(loop);
   lamda6=lamda2(loop);
  
% COMPUTATION OF DISPLACEMENT SENSITIVITIES
 [U5]=reddistrib(nelx,nely,x,dU,dU1,lamda5,lamda6,p,a,b,h);
 % COMPUTATION OF FREQUENCE SENSITIVITIES
 [G,G1,G2]=refdistrib(nelx,nely,G,p,a,b,dU,dU1);
% DESIGN UPDATE BY THE OPTIMALITY CRITERIA METHOD
 [x]=OC(nelx,nely,x,h,U5,q,G,freq,d,g,freq1,p,a,b,dU,dU1,G1,G2);  
 % m=1-(1-0.05)^((loop+1)/loop);
%  m=0.2+0.5^((loop+1)/loop);
y(loop)=c;
% PRINT RESULTS
change = abs(4*p*a*b*sum(sum(x))-4*p*a*b*sum(sum(xold)))/abs(4*p*a*b*sum(sum(x)));
  disp([' It.: ' sprintf('%4i',loop) ' Obj.: ' sprintf('%10.4f',4*p*a*b*sum(sum(x))) ...
       ' Vol.: ' sprintf('%6.3f',4*sum(sum(x))*a*b) ...
        ' lamda5.: ' sprintf('%6.3f',lamda5 )])
% PLOT DENSITIES  
figure(1)
  colormap(gray); imagesc(-x); axis equal; axis tight; axis off;pause(1e-6);
figure(2)
plot(y)
end 
