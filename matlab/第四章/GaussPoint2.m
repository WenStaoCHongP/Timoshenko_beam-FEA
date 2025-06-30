%--------------------------------------------------------------------------
function [point2,weight2]=GaussPoint2(No_point_x,No_point_y)
%--------------------------------------------------------------------------
%  Purpose:
%     determine the integration points and weighting coefficients
%     of Gauss quadrature for two-dimensional integration
%--------------------------------------------------------------------------
%----------------------------------------------------------
%  (1) determine the largest one between No_point_x and No_point_y
%----------------------------------------------------------
   if No_point_x > No_point_y
      ng=No_point_x;
   else
      ng=No_point_y;
   end
%----------------------------------
%  (2) initialization of coefficient vector
%----------------------------------
   point2=zeros(ng,2);
   weight2=zeros(ng,2);
%-----------------------------------------------
%  (3) find corresponding integration points and weights
%-----------------------------------------------
 [pointx,weightx]=GaussPoint1(No_point_x);     % quadrature rule for x-axis
 [pointy,weighty]=GaussPoint1(No_point_y);     % quadrature rule for y-axis
%-------------------------------
%  ги4) quadrature for two-dimension
%-------------------------------
 for ii=1:No_point_x                         % quadrature in x-axis
   point2(ii,1)=pointx(ii);
   weight2(ii,1)=weightx(ii);
 end
 for ij=1:No_point_y                         % quadrature in y-axis
   point2(ij,2)=pointy(ij);
   weight2(ij,2)=weighty(ij);
 end
%--------------------------------------------------------------------------
%    The end
%--------------------------------------------------------------------------
