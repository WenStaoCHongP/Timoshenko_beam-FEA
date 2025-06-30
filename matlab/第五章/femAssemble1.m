%--------------------------------------------------------------------------
function [kk]=femAssemble1(kk,k,index)
%--------------------------------------------------------------------------
%  Purpose:
%     Assembly of element matrices into the system matrix
%  Synopsis:
%     [kk]=femAssemble1(kk,k,index)
%  Variable Description:
%     kk - system matrix
%     k  - element matrix
%     index - dof vector associated with an element
%--------------------------------------------------------------------------
eldof = length(index);

for i=1:eldof
   ii=index(i);
   for j=1:eldof
      jj=index(j);
      kk(ii,jj)=kk(ii,jj)+k(i,j);
   end
end
%--------------------------------------------------------------------------
%     The end
%--------------------------------------------------------------------------
