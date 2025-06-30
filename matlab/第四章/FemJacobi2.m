%--------------------------------------------------------------------------
function [jacob2]=FemJacobi2(No_nodeEl,dhdr,dhds,xcoord,ycoord)
%--------------------------------------------------------------------------
%  Purpose:
%     determine the Jacobian for two-dimensional mapping
%  Variable Description:
%     jacob2 - Jacobian for one-dimension
%     No_nodeEl - number of nodes per element
%     dhdr - derivative of shape functions w.r.t. natural coordinate r
%     dhds - derivative of shape functions w.r.t. natural coordinate s
%     xcoord - x axis coordinate values of nodes
%     ycoord - y axis coordinate values of nodes
%------------------------------------------------------------------------
 jacob2=zeros(2,2);
 
 for i=1: No_nodeEl
 jacob2(1,1)=jacob2(1,1)+dhdr(i)*xcoord(i);
 jacob2(1,2)=jacob2(1,2)+dhdr(i)*ycoord(i);
 jacob2(2,1)=jacob2(2,1)+dhds(i)*xcoord(i);
 jacob2(2,2)=jacob2(2,2)+dhds(i)*ycoord(i);
 end
%--------------------------------------------------------------------------
%    The end
%--------------------------------------------------------------------------
