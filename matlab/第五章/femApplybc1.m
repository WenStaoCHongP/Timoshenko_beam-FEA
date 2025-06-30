%--------------------------------------------------------------------------
function [kk,mm,ff]=femApplybc1(kk,mm,ff,bcdof,bcval)
%--------------------------------------------------------------------------
%  Purpose:
%     Apply constraints to matrix equation [kk]{x}={ff}
%     Apply constraints to eigenvalue matrix equation 
%     [kk]{x}=lamda[mm]{x}
%  Synopsis:
%     [kk,mm,ff]=femApplybc1(kk,mm,ff,bcdof,bcval)
%  Variable Description:
%     kk - system stiffness matrix before applying constraints 
%     mm - system mass matrix before applying constraints
%     ff - system vector before applying constraints
%     bcdof - a vector containing constrained dof
%     bcval - a vector containing contained value 
%--------------------------------------------------------------------------
 ni=length(bcdof);
 sdof=size(kk);
 
 for ii=1:ni
   if bcdof(ii)==1
     for j=1:sdof
       kk(ii,j)=0;
       kk(j,ii)=0;
       mm(ii,j)=0;
       mm(j,ii)=0;
       ff(j)=ff(j)-bcval(ii)*kk(j,ii);
     end
     kk(ii,ii)=1;
     ff(ii)=bcval(ii);
   end
 end
%--------------------------------------------------------------------------
%     The end
%--------------------------------------------------------------------------
