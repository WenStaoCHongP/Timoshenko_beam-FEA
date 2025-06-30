%-------------------------------------------------------------------------------------
 function [k, m]=fetruss31(leng,c,s,r)

%-------------------------------------------------------------------------------------
%  Purpose:
%     Stiffness and mass matrices for piezpelectric element 
%
%  Synopsis:
%     [k,m]=fetruss31(leng,c,s,r) 
%
%  Variable Description:
%     k - element stiffness matrix (size of 6x6)   
%     m - element mass matrix (size of 6x6)
%    leng - element length
%     c,s,r -  direction cosine between global coordinate and local coordinate
%    
%   One assumption: the lengths of two connected bars are same
%   
%--------------------------------------------------------------------------
 el1=210e9;                    % Young's module of connected bars
 rho1=8000;                    % mass density of connected bars
 area1=6.9e-5;                  % cross-section area of connected bars
 el2=63e9;                     % Young's module of piezoelectric stack
 rho2=7600;                    % mass density of piezoelectric stack
 area2=7.07e-4;                 % cross-section area of piezoelectric stack
 leng2=0.56;
 k1=area1*el1/((leng-leng2)/2); 
 k2=(area2*el2)/leng2;
 k12=(k1^2*k2)/(2*k1*k2+k1^2);     
 m12=rho1*area1*(leng-leng2)+rho2*area2*leng2;   % total mass of piezoelectric active member
% stiffness matrix 

 k13= k12*[ c*c c*s c*r;...
            c*s s*s s*r;...
           c*r s*r r*r];

 k=[k13 -k13;-k13 k13];
% lumped mass matrix

 m=m12/3*eye(6);
%
%-------------------------------------------------------------------------------------
