%-------------------------------------------------------------
% Ex.8.1.5
%  one side with fixly-supported condition for rectangular plate with
%  individual displacement and frequence constraint
%  30*16 Mesh
%-------------------------------------------------------------
% INITIALIZE
nelx=52;
nely=26;
x(1:nely,1:nelx) = 0.005*ones(nely,nelx); 
a=0.005;
b=0.005;
p=1e6;
loop = 0; 
h=0.21;
U44=-0.0005;
freq1=70;
change =0.01;
% START ITERATION
q=0.2;
d=1.4;       % parameter 2
g=0.28;       % parameter 3
while change>0.00001
 loop = loop + 1;
  xold = x; 
% FE-ANALYSIS
 
  [U,U2,G,freq,lambda]=FE(nelx,nely,xold); 
  if ((loop)>120),break,end
% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
 
[KE,ME] = lk;
    c = 0.;
    U99=U(2*(0.5*nelx+1)*(nely+1),1)
    freq(1)
  for ely = 1:nely
    for elx = 1:nelx
      n1 = (nely+1)*(elx-1)+ely; 
      n2 = (nely+1)* elx   +ely;
      edof = [2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n2+1; 2*n2+2; 2*n1+1; 2*n1+2];
      c = c + 4*((xold(ely,elx)))*p*a*b;
      dc(ely,elx) = 4*p*a*b;
      Ue = U([2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n2+1; 2*n2+2; 2*n1+1; 2*n1+2],1);
      F2(edof,:)=xold(ely,elx)*KE*Ue;
      dU(ely,elx)=F2(edof,:)'*U2([edof],2*(0.5*nelx+1)*(nely+1));
      end  
  end
% Lagrange multiplier
    lamda=c/(3*abs(U44));
% COMPUTATION OF DISPLACEMENT SENSITIVITIES
 [U5]=reddistrib(nelx,nely,x,dU,lamda,p,a,b,h);
 [G,G1,G2]=refdistrib(nelx,nely,G,p,a,b,dU);
% DESIGN UPDATE BY THE OPTIMALITY CRITERIA METHOD
[x]=OC(nelx,nely,x,h,U5,q,G,freq,d,g,freq1,p,a,b,dU,G1,G2); 

% iterated objective function
y(loop)=c;

% PRINT RESULTS
change = abs(4*p*a*b*sum(sum(x))-4*p*a*b*sum(sum(xold)))/abs(4*p*a*b*sum(sum(x)));
  disp([' It.: ' sprintf('%4i',loop) ' Obj.: ' sprintf('%6.4f',4*p*a*b*sum(sum(x))) ...
       ' Vol.: ' sprintf('%6.3f',4*sum(sum(x))*a*b) ...
        ' lamda.: ' sprintf('%6.3f',lamda )])

% PLOT DENSITIES

figure(1)
  colormap(gray); imagesc(-x); axis equal; axis tight; axis off;pause(1e-6); 
figure(2)
plot(y)
end
