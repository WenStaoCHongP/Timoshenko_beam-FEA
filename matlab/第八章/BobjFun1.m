function [f]=BobjFun1(x)
%--------------------------------------------------------------------------
%  Purpose:                                                               
%     Calculate the objective function which is the norm of difference of 
%     modal vectors                                                                                                                             
%  Synopsis:                                                              
%     [f]=BobjFun1(x)                                                                                                                          
%  Variable Description:                                                  
%     x = design variables                                                
%     f = objective function                                              
%-------------------------------------------------------------------------------
%  (1) prescribe frequency and ith order modal verctor 
%-------------------------------------------------------------------------------
num_f=5;                        %  number of natural frequency
ith_m=2;                        %  1-ith order of modal vectors
kg=0;                           % option for graphic display of modal shape
global Vit
%------------------------------------------------------------------------------
%  (2) calculte ith order modal verctor of two structure
%------------------------------------------------------------------------------
[freq,ws,Vi]=BeamCom1(x,num_f,ith_m,kg);                  % FEM calculation
                            % freq == natrural requency of structure system
                            % ws == total weight of boundary elements
                            % Vi == modal vectors
%-------------------------------------------------------------------------------------
%  (3) calculate the norm of difference between two modal vectors
%-------------------------------------------------------------------------------------
f=sum((Vi-Vit).*(Vi-Vit));
f=sum(f);
f=sqrt(f);
% -------------------------------------------------------------------------
%     The end
% -------------------------------------------------------------------------
