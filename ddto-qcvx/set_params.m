%{
06/04/2022
Purnanand Elango

DDTO algorithm for 3D double integrator with gravity
%}

global nx nu Ad Bd u_max u_min stage_cost d alf_tp

Del = 0.5; % Sampling time

cvx_solver mosek

nx = 6; % System dimension
nu = 3; % Number of control inputs

m = 1;
g = 9.806;
alf_tp = tand(60);

z0 = [0; 0; 30; zeros(3,1)]; % Initial state

% Target states
zf = 40*[real(exp(1i*pi/20))   real(exp(1i*pi/4))  real(exp(1i*pi/2))  real(exp(-1i*pi/20));
          imag(exp(1i*pi/20))   imag(exp(1i*pi/4))  imag(exp(1i*pi/2))  imag(exp(-1i*pi/20));
          zeros(4,4)];

N = [40,40,40,40]; % Horizon length
  
prob_par = struct;  
prob_par.z0 = z0;                   % Initial state
prob_par.zf = zf;                   % Target states   
prob_par.N = N;                     % No. of discretization points 
prob_par.n = size(zf,2);            % No. of targets   
prob_par.tag = 1:length(zf(1,:));   % Target tag
prob_par.opt_tol = 0.2;             % Optimality tolerance

clear z0 zf N

u_min = 8; % Input z-axis component lowerbound
u_max = 20; % Input norm upperbound

% Define system matrices
Ad = [1     0     0     Del     0      0;
      0     1     0      0     Del     0;
      0     0     1      0      0     Del;
      0     0     0      1      0      0;
      0     0     0      0      1      0;
      0     0     0      0      0      1];

Bd = [Del*Del/(2*m)        0             0;
           0          Del*Del/(2*m)      0;
           0               0        Del*Del/(2*m);   
         Del/m             0             0;
           0             Del/m           0;
           0               0           Del/m];

d = [0; 0; -0.5*Del*Del*g; 0; 0; -Del*g];
     
% Cost weights
Q = eye(nx);
R = eye(nu);
     
% stage_cost = @(x,u) norm(u,2);
stage_cost = @(x,u) u'*R*u;
