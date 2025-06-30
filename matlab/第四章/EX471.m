%--------------------------------------------------------------------------
% Example 4.7.1
%--------------------------------------------------------------------------
%   Gauss quadrature of a function in one, two, three dimension
% Problem descriptions
%   Integrate:
%   (1)  f(x)=1+x^2-3*x^3+4*x^5  (-1<x<1)
%   (2)  f(x,y)=1+4*x*y-3*x^2*y^2+x^4*y^6   (-1<x<1, -1<y<1)
%   (3)  f(x,y,z)=1+4*x^2*y^2-3*x^2*z^4+y^4*z^6  (-1<(x,y,z)<1)
% Variable descriptions
%   points = integration (or sampling) points
%   weights = weighting coefficients
%   No_INTpoint_x = number of integration points along x-axis
%   No_INTpoint_t = number of integration points along y-axis
%   No_INTpoint_z = number of integration points along z-axis
%--------------------------------------------------------------------------
%----------------------------
%  (1) Selecte the problem
%----------------------------
clear; clc;

Opt_problem=1;                           % option for type of the beam:
                                        % =1 One-dimensional function
                                        % =2 Two-dimensional function
                                        % =3 Three-dimensional function
No_INTpoint_x=5;                         % No_INTpoint_x >= 0.5*(Order_x+1)
No_INTpoint_y=0;                         % No_INTpoint_y >= 0.5*(Order_y+1)
No_INTpoint_z=0;                         % No_INTpoint_z >= 0.5*(Order_z+1)
%-------------------------------
%  (2) initialize the vectors
%-------------------------------
   pointx=zeros(No_INTpoint_x,1);
   weightx=zeros(No_INTpoint_x,1);

   pointy=zeros(No_INTpoint_y,1);
   weighty=zeros(No_INTpoint_y,1);

   pointz=zeros(No_INTpoint_z,1);
   weightz=zeros(No_INTpoint_z,1);
%--------------------------------------------
%  (3) compute numerical integration
%--------------------------------------------
value=0.0;

switch Opt_problem
  case 1
    [pointx,weightx]=GaussPoint1(No_INTpoint_x);
                                      % extract integration points and weights
    for ii=1:No_INTpoint_x
      x=pointx(ii);                      % sampling points in x-axis
      wtx=weightx(ii);                   % weighting coefficients in x-axis
      func=1+x^2-3*x^3+4*x^5;                           % evaluate function
      value=value+func*wtx;
    end
  case 2
    [pointx,weightx]=GaussPoint1(No_INTpoint_x);
    [pointy,weighty]=GaussPoint1(No_INTpoint_y);
                                      % extract integration points and weights
    for ii=1:No_INTpoint_x
      x=pointx(ii);                      % sampling points in x-axis
      wtx=weightx(ii);                   % weighting coefficients in x-axis
      for ij=1:No_INTpoint_y
        y=pointy(ij);                    % sampling points in y-axis
        wty=weighty(ij) ;                % weighting coefficients in y-axis
        func=1+4*x*y-3*x^2*y^2+x^4*y^6;                 % evaluate function
        value=value+func*wtx*wty;
      end
    end
  case 3
    [pointx,weightx]=GaussPoint1(No_INTpoint_x);
    [pointy,weighty]=GaussPoint1(No_INTpoint_y);
    [pointz,weightz]=GaussPoint1(No_INTpoint_z);
                                      % extract integration points and weights
    for ii=1:No_INTpoint_x
      x=pointx(ii);                      % sampling points in x-axis
      wtx=weightx(ii);                   % weighting coefficients in x-axis
      for ij=1:No_INTpoint_y
        y=pointy(ij);                    % sampling points in y-axis
        wty=weighty(ij) ;                % weighting coefficients in y-axis
        for ik=1:No_INTpoint_z
          z=pointz(ik);                  % sampling points in z-axis
          wtz=weightz(ik) ;              % weighting coefficients in z-axis
          func=1+4*x^2*y^2-3*x^2*z^4+y^4*z^6;           % evaluate function
          value=value+func*wtx*wty*wtz;
        end
      end
    end
  otherwise
    disp('wrong number for problem')
end
 
value                                  % print the solution
%--------------------------------------------------------------------------
%    The end
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
function [point1,weight1]=GaussPoint1(No_points)
%--------------------------------------------------------------------------
%  Purpose:
%     determine the integration points and weighting coefficients of Gauss quadrature
%  Variable Description:
%     No_points - number of integration points
%     point1  - vector containing integration points
%     weight1 - vector containing weighting coefficients
%-------------------------------------------------------------------
%------------------------------------------
%  (1) initialization of the vectors
%------------------------------------------
   points=zeros(No_points,1);
   weights=zeros(No_points,1);
%----------------------------------------------------------
%  (2) find corresponding integration points and weights
%----------------------------------------------------------
if No_points==1                 % 1-point quadrature rule
    point1(1)=0.0;
    weight1(1)=2.0;
elseif No_points==2             % 2-point quadrature rule
    point1(1)=-0.577350269189626;
    point1(2)=-point1(1);
    
    weight1(1)=1.0;
    weight1(2)=weight1(1);
elseif No_points==3             % 3-point quadrature rule
    point1(1)=-0.774596669241483;
    point1(2)=0.0;
    point1(3)=-point1(1);
    
    weight1(1)=0.555555555555556;
    weight1(2)=0.888888888888889;
    weight1(3)=weight1(1);
elseif No_points==4             % 4-point quadrature rule
    point1(1)=-0.861136311594053;
    point1(2)=-0.339981043584856;
    point1(3)=-point1(2);
    point1ги4)=-point1(1);
    
    weight1(1)=0.347854845137454;
    weight1(2)=0.652145154862546;
    weight1(3)=weight1(2);
    weight1ги4)=weight1(1);
elseif No_points==5             % 5-point quadrature rule
    point1(1)=-0.906179845938664;
    point1(2)=-0.538469310105683;
    point1(3)=0.0;
    point1ги4)=-point1(2);
    point1(5)=-point1(1);
    
    weight1(1)=0.236926885056189;
    weight1(2)=0.478628670499366;
    weight1(3)=0.568888888888889;
    weight1ги4)=weight1(2);
    weights(5)=weight1(1);
elseif  No_points==6            % 6-point quadrature rule
    point1(1)=-0.932469514203152;
    point1(2)=-0.661209386466265;
    point1(3)=-0.238619186083197;
    point1ги4)=-point1(3);
    point1(5)=-point1(2);
    point1(6)=-point1(1);
    
    weight1(1)=0.171324492379170;
    weight1(2)=0.360761573048139;
    weight1(3)=0.467913934572691;
    weight1ги4)=weight1(3);
    weight1(5)=weight1(2);
    weight1(6)=weight1(1);
else
   disp('number of integration points should be')
   disp('> 0 or < 7')
end
%-------------------------------------------------------------------
%    The end
%-------------------------------------------------------------------
