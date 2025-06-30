%--------------------------------------------------------------------------
function [index]=femEldof(nd,No_nel,No_dof)
%--------------------------------------------------------------------------
%  Purpose:
%     Compute system dofs associated with each element
%  Synopsis:
%     [index]=femEldof(nd,No_nel,No_dof)
%  Variable Description:
%     index - system dof vector associated with element "iel"
%     nd - node connectivity for the (iel)-th element
%     No_nel - number of nodes per element
%     No_dof - number of dofs per node
%     iel - element number whose system dofs are to be determined
%--------------------------------------------------------------------------
   k=0;
   for i=1:No_nel
     start = (nd(i)-1)*No_dof;
       for j=1:No_dof
         k=k+1;
         index(k)=start+j;
       end
   end
%--------------------------------------------------------------------------
%     The end
%--------------------------------------------------------------------------
