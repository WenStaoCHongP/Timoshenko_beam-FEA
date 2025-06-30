%-------------------------------------------------------------------------------------
% BeamProbB21  (minimize the difference of lower order modal vectors)     
%   BeamProb1: (4 supports)                                               
%   To find the sizes of the structure using optimal algorithm            
%                                                                         
% Variable descriptions                                                   
%   Inout:                                                                
%         x0 == intial value of the design variables                      
%   Optput:                                                               
%         x == desin varialbes at the optimal solution                    
%         fval == value of objective function                             
%--------------------------------------------------------------------------
%  (0) input control data
%------------------------------------
clear; clc;
num_f=5;                        %  number of natural frequency
ith_m=2;                         %  1-ith order of modal vectors
kg=0;                           % option for graphic display of modal shape
global freqt wst Vit
                  % freqt == natrural requency of oririnal structure system
                  % wst == total weight of boundary elements
                  % Vit == modal vector of oririnal structure system
%----------------------------------------------------------------------------------------
%  (1) calculate the dynamic properties of original structure
%--------------------------------------------------------------------------------------
x0=[0.0148; 0.0148; 0.0148; 0.0148; 0.0148; 0.0148;
    0.0040; 0.0040; 0.0040; 0.0040; 0.0040; 0.0040];
[freqt,wst,Vit]=BeamCom0(x0,num_f,ith_m,kg);              % FEM calculation
%-----------------------------------------
%  (2) input variable data
%------------------------------------------
                               % intial value of number of design variables
h0=0.0148; b0=0.005;
x0=[h0; h0; h0; h0;
    b0; b0; b0; b0];              % set the lower and upper bounds of the variables
%hl=0.004; hu=0.020;
%bl=0.004; bu=0.020;
hl=0.0148; hu=0.0148;
bl=0.0040; bu=0.0100;
vlb=[hl; hl; hl; hl;
     bl; bl; bl; bl];
vub=[hu; hu; hu; hu;
     bu; bu; bu; bu];
%--------------------------------------------------------------------------
%  (3) process the optimum numerical operation
%-------------------------------------------------------------------------
% call the fmincon function to get optimal solution                  
tic
options = optimset('LargeScale','off');
options = optimset('MaxSQPIter', 1000);
options = optimset('MaxFunEvals',2000);
[x,fval]=fmincon(@BobjFun1,x0,[],[],[],[],vlb,vub,@ConFun1,options)
[c, ceq, freq]=ConFun1(x)
te=toc
%--------------------------------------------------------------------------
%     The end
% -------------------------------------------------------------------------
