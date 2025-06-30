%-------------------------------------------------------------------------%
% Example 5.4.8
%   To solve transient response of a 2-d truss structure and
%   calculate Fast Fourier Transform (FFT) of the time domain response.
%   The solution methods are: 1) central difference scheme. 3) Houbolt integration scheme.
%   4) Wilson   integration scheme. 5) Newmark integration scheme
%   nodal dof: {u1  v1  w1  x1  y1  z1  u2  v2  w2  x2  y2  z2}
% Problem description
%   Find the response of a frame structure which is made of three beams of lengths of 4 m, 
%   3 m and 4 m, respectively. All beams have cross- section of 0.10 m height by 0.05 m
%   width. The elastic modulus is 2.10x10^11 Pa. The frame is subjected to an impulse load
%   of amplitude 500 N in the middle of the top beam. One end of the each vertical beam is
%   fixed. (see Fig. 5-9 for the element discretization)
% Variable descriptions
%   k, m - element stiffness matrix and mass matrix
%   kk, mm - system stiffness matrix and mass matrix
%   ff - system force vector
%   index - a vector containing system dofs associated with each element
%   bcdof - a vector containing dofs associated with boundary conditions
%   bcval - a vector containing boundary condition values associated with the dofs in 'bcdof'
%   dsp - displacement matrix
%   vel - velocity matrix
%   acc - acceleartion matrix
%--------------------------------------------------------------------------
%  (0) input control data
%--------------------------------------------------------------------------
clear; clc;

Beam_InputData547;                     % import the input data for the information of
                                  % nodes, elements, loads, constraints and materials
Opt_beam=1;                                        % option for type of the beam
                                                     % =1 Euler Bernoulli beam
                                                       % =2 Timoshenko beam
Opt_mass=2;                                            % option for mass matrix
                                                    % =1 consistent mass matrix
                                                      % =2 lumped mass matrix
Opt_section=1;                                    % option for type of cross-section
                                                  % = 1 rectangular cross-section
                                                     % = 2 circular cross-section
TypeMethod=4;                            % option for selecting the solution method
                                                 % = 1 central difference scheme
                                               % = 3 Houbolt integration scheme
                                             % = 4 Wilson   integration scheme
                                              % = 5 Newmark integration scheme
Typeload=1;                                    % option for selecting the load type
                                                           % = 1 impulse load
                                                              % = 2 step load
                                                          % = 3 Harmonic load

dt=0.0005;                                                     % time step size
ti=0;                                                            % initial time
tf=0.200;                                                          % final time
nt=fix((tf-ti)/dt);                                           % number of time steps
tt=ti:dt:ti+nt*dt;                                     % generate time samples vector
 
ac=0.00002; bc=0.00008;                       % Parameters for proportional damping
 
al=0;                            % angle between the reference coordinate system and
                                  % the local coordinate system for the space element
%--------------------------------------------------------------------------
%  (7) calculate transient response
%--------------------------------------------------------------------------
[kk1,mm1,ff1,bcdof1,bcval1,sdof1]=mmCheck1(kk0,mm0,ff0,bcdof0,bcval0);
                         % check the zero main elements in mm and eliminate the rows
                      % and columns in equation associated with the zero main elements
switch Typeload
  case 1                                                % Impulse force function
    u=[1,zeros(1,nt)];
    ft0=ff1*u;
  case 2                                                   % Step force function
    u(1,1:nt+1)=1;
    ft0=ff1*u;
  case 3                                              % Harmonic force function
    u=cos(omega0*tt);
    ft0=ff1*u;
  otherwise
    ft0=ff1;                                            % a given force function
end

cc1=ac*mm1+bc*kk1;                          % Form the Rayleigh damping matrix
q0=zeros(sdof1,1); dq0=zeros(sdof1,1);               % initial displacement and velocity

switch TypeMethod
  case 1                                              % central difference scheme
    [acc,vel,dsp]=TransResp1(kk1,cc1,mm1,ft0,bcdof1,nt,dt,q0,dq0);
  case 3                                            % Houbolt integration scheme
    [acc,vel,dsp]=TransResp3(kk1,cc1,mm1,ft0,bcdof1,nt,dt,q0,dq0);
  case 4                                         % Wilson   integration scheme
    [acc,vel,dsp]=TransResp4(kk1,cc1,mm1,ft0,bcdof1,nt,dt,q0,dq0);
  case 5                                           % Newmark integration scheme
    [acc,vel,dsp]=TransResp5(kk1,cc1,mm1,ft0,bcdof1,nt,dt,q0,dq0);
  otherwise
    disp('Unknown method.')
end
%--------------------------------------------------------------------------
%  (8) Calculate FFT of the time domain data y(n,:)
%--------------------------------------------------------------------------
jth=20;
 
yt=dsp(jth,:)';                        % extract one time domain data from the response
tt=ti:dt:ti+nt*dt;                                            % generate time vector
 
[yfft, freq]=femFFT(yt,tt);                %  Calculate FFT of the time domain data and
                                              %  take absolute values of the result
%--------------------------------------------------------------------------
%  (9) graphics of dynamic response
%--------------------------------------------------------------------------
jth=20;

figure(1)
                    % Plot the graph of the time-history response at jth degree of freedom
plot(tt,dsp(jth,:))
xlabel('Time  (seconds)'), ylabel('displacement  (m)')
title('time-history response')
 
figure(2)                                % Plot the graph of the FFT versus frequency
plot(freq,yfft)
xlabel('Frequency (rad/s)')
ylabel('Absolute values of FFT')
title('FFT result')
%--------------------------------------------------------------------------
%  (10) print fem solutions
%--------------------------------------------------------------------------
switch Typeload
  case 1
    disp('The excitation is impulse force')
  case 2
    disp('The excitation is step force')
  case 3
    disp('The excitation is step force')
    otherwise
    disp('The given foece')
end

disp('The calculation is use of:')

if Opt_beam==1
  disp('Euler-Bernoulli beam element 1')
elseif Opt_beam==2
  disp('Timoshenko beam element')
else
  disp('Euler-Bernoulli beam element 2')
end
 
if Opt_mass==1
  disp('and consistent mass matrix')
elseif Opt_mass==2
  disp('and lumped mass matrix')
else
  disp('and diagonal mass matrix')
end
 
num=1:1:sdof2;
frequency=[num' omega1]                                  % print natural frequency
%--------------------------------------------------------------------------
%     The end
%--------------------------------------------------------------------------
