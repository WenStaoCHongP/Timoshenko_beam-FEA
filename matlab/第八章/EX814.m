%-------------------------------------------------------------------------------------
% Ex.8.1.4
% one side with fixly-supported condition for rectangular plate with individual 
%displacement and  frequency constraints
%-------------------------------------------------------------------------------------
% INITIALIZE
nelx=20;   % the numbers of elements of x direction
nely=48;   % the numbers of  elements of y direction
x(1:nely,1:nelx) = 0.006*ones(nely,nelx);       %  initialization of thickness of elements
a=0.0025;  % length of element along x direction
b=0.0025;  % length of element along y direction
p=1e4;     %  density of material
loop = 0;   % iteration variable
h=0.2;     % relaxation factor
U44=-0.00015;  % constrained displacement 
freq1=3000;    % constrained frequency 
change =0.01;  
q=0.21;      % parameter1
d=1.4;       % parameter2 about frequency sensitivity
g=0.1;       % parameter3 about frequency sensitivity
% START ITERATION
while change>0.0001
 loop = loop + 1;
  xold = x; 
% FE-ANALYSIS
  [U,U2,G,freq,lambda]=FE(nelx,nely,xold); 
% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
 
[KE,ME] = lk;
    c = 0.;
  for ely = 1:nely
    for elx = 1:nelx
      n1 = (nely+1)*(elx-1)+ely; 
      n2 = (nely+1)* elx   +ely;
      edof = [2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n2+1; 2*n2+2; 2*n1+1; 2*n1+2];
      c = c + 4*((xold(ely,elx)))*p*a*b;
      dc(ely,elx) = 4*p*a*b;
      Ue = U([2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n2+1; 2*n2+2; 2*n1+1; 2*n1+2],1);
      F2(edof,:)=xold(ely,elx)*KE*Ue;
      dU(ely,elx)=F2(edof,:)'*U2([edof],2*(nelx+0.5)*(nely+1)+1);
      end  
  end
% Lagrange multiplier
    lamda=c/(3*abs(U44));
% COMPUTATION OF DISPLACEMENT SENSITIVITIES
 [U5]=reddistrib(nelx,nely,x,dU,lamda,p,a,b,h);
% COMPUTATION OF FREQUENCY SENSITIVITIES
 [G,G1,G2]=refdistrib(nelx,nely,G,p,a,b,dU);
% DESIGN UPDATE BY THE OPTIMALITY CRITERIA METHOD
[x]=OC(nelx,nely,x,h,U5,q,G,freq,d,g,freq1,p,a,b,dU,G1,G2); 
% iterated objective function
y(loop)=c;
% PRINT RESULTS
change = abs(4*p*a*b*sum(sum(x))-4*p*a*b*sum(sum(xold)))/abs(4*p*a*b*sum(sum(x)));
  disp([' It.: ' sprintf('%4i',loop) ' Obj.: ' sprintf('%10.4f',4*p*a*b*sum(sum(x))) ...
       ' Vol.: ' sprintf('%6.3f',4*sum(sum(x))*a*b) ...
        ' freq(1).: ' sprintf('%6.3f',freq(1) )])

% PLOT DENSITIES
figure(1)
  colormap(gray); imagesc(-x); axis equal; axis tight; axis off;pause(1e-6); grid on;
figure(2)
plot(y)
end
