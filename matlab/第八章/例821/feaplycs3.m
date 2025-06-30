function [kkm, mmm]=feaplycs3(kk, mm, bcdof)
%-------------------------------------------------------------------------------------
% Purpose:
% Apply constraints to obtain the final system mass matrix and stiffness matrix 
%
%  Synopsis:
%     [kkm,mm]=feaplycs(kk,mm,bcdof)
%
%  Variable Description:
%    kkm- system stiffness matrix after applying constraints
%    mmm- system mass matrix after applying constraints
%     kk - system stiffness matrix before applying constraints 
%     mm - system mass matrix before applying constraints
%     bcdof - a vector containging constrained d.o.f
%-------------------------------------------------------------------------------------bcdof=sort(bcdof);
n=length(bcdof);
j=-1
for i=1:n
j=j+1;
mm(:, (bcdof(i)-j))=[];
mm((bcdof(i)-j), :)=[];
kk(:, (bcdof(i)-j))=[];
kk((bcdof(i)-j), :)=[];
end
kkm=kk;
mmm=mm;
%--------------------------------------------------------------------
