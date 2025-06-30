function [c, ceq, freq] = ConFun1(x)
%--------------------------------------------------------------------------
%  Purpose:                                                               
%     To computes the nonlinear inequality constraints c(x)<= 0 and the   
%     nonlinear equality constraints ceq(x) = 0. This function accepts    
%     a vector x and returns two vectors c and ceq. The vector c contains 
%     the nonlinear inequalities evaluated at x, and ceq contains the     
%     nonlinear equalities evaluated at x.                                                                                                        
%  Synopsis:                                                              
%     [c, ceq]=ConFun1(x)                                                                                                                     
%  Variable Description:                                                  
%    Input:                                                               
%          x == design variables                                          
%    Optput:                                                              
%          c == inequality constraint conditions                          
%          ceq == equality constraint conditions                          %
%--------------------------------------------------------------------------
%  (0) input control data
%----------------------------------------
num_f=5;                       %  number of natural frequency
ith_m=2;                        %  1-ith order of modal vectors
kg=0;                           % option for graphic display of modal shape
global freqt
%-------------------------------------------------------------------------------------
%  (1) calculte the frequency and 1-ith order modal verctors
%-------------------------------------------------------------------------------------
[freq,ws,Vi]=BeamCom1(x,num_f,ith_m,kg);                  % FEM calculation
                            % freq == natrural requency of structure system
                            % ws == total weight of boundary elements
                            % Vi == modal vectors
%--------------------------------------------------------------------------
%  (2) check the boundary conditions
%--------------------------------------------------------------------------
                                         % Nonlinear inequality constraints
for i=1:num_f
  c(i)=abs(freq(i)-freqt(i))-0.01*freqt(i);
end    
                                        % No nonlinear equality constraints
ceq=[];
% -------------------------------------------------------------------------
%     The end
% -------------------------------------------------------------------------
