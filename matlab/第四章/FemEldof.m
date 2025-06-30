%--------------------------------------------------------------------------
function [index]=FemEldof(nd, No_nodeEl,No_dof)
%--------------------------------------------------------------------------
%  Purpose:
%     Compute system dofs associated with each element
%  Variable Description:
%     index - system dof vector associated with element "iel"
%     nd - node connectivity for the (iel)-th element
%     No_nodeEl - number of nodes per element
%     No_dof - number of dofs per node
%--------------------------------------------------------------------------
   k=0;
   for i=1: No_nodeEl
     start = (nd(i)-1)* No_dof;
       for j=1:No_dof
         k=k+1;
         index(k)=start+j;
       end
   end
%--------------------------------------------------------------------------
%    The end
%--------------------------------------------------------------------------
