function prb = problem_data(preferred_target,M)

prb = struct;

prb.M = M;                                                                                                      % Constant in big-M reformulation of MIP

prb.Del = 0.5;                                                                                                  % Sampling time [s]

prb.nx = 6;                                                                                                     % State dimension
prb.nu = 3;                                                                                                     % Control input dimension

prb.mass = 1;                                                                                                   % Vehicle mass [kg]
prb.accl = 9.806;                                                                                               % Acceleration due to gravity [m/s^2]
prb.tan_thet_tp = tand(60);                                                                                     % Tangent of thrust pointing angle     

prb.z0 = [0; 0; 30; zeros(3,1)];                                                                                % Initial state [m,m/s]

                                                                                                                % Target states [m,m/s]
prb.zf = 40*[real(exp(1i*pi/20))   real(exp(1i*pi/4))  real(exp(1i*pi/2))  real(exp(-1i*pi/20));
             imag(exp(1i*pi/20))   imag(exp(1i*pi/4))  imag(exp(1i*pi/2))  imag(exp(-1i*pi/20));
             zeros(4,4)];

prb.N = 20*[1,1,1,1];                                                                                           % Trajectory length
% prb.N = [20,25,22,23];

prb.n = size(prb.zf,2);                                                                                         % No. of targets

prb.i = preferred_target;                                                                                       % Preferred target

prb.Ni = zeros(1,prb.n-1);                                                                                      % min{Nj,Ni}
for j = setdiff(1:prb.n,prb.i)
    prb.Ni(j) = min(prb.N(j),prb.N(prb.i));
end

prb.lmax = 3800; %2150                                                                                                % Cumulative trajectory cost upper bound    

prb.umin = 8;                                                                                                   % Lower bound on control input z-axis component
prb.umax = 20;                                                                                                  % Upper bound on control input magnitude

% Define system matrices
prb.A  = [1     0     0      prb.Del     0          0;
          0     1     0      0           prb.Del    0;
          0     0     1      0           0          prb.Del;
          0     0     0      1           0          0;
          0     0     0      0           1          0;
          0     0     0      0           0          1];

prb.B  = [prb.Del*prb.Del/(2*prb.mass)  0                               0;
          0                             prb.Del*prb.Del/(2*prb.mass)    0;
          0                             0                               prb.Del*prb.Del/(2*prb.mass);   
          prb.Del/prb.mass              0                               0;
          0                             prb.Del/prb.mass                0;
          0                             0                               prb.Del/prb.mass];

prb.c  = [ 0; 
           0; 
          -0.5*prb.Del*prb.Del*prb.accl; 
           0; 
           0; 
          -prb.Del*prb.accl];
     
prb.stage_cost = @(x,u) u'*u;                                                                                   % Stage cost for defining cumulative trajectory cost

% prb.solversettings = sdpsettings('solver','gurobi','verbose',0,'gurobi.OptimalityTol',1e-9,'gurobi.FeasibilityTol',1e-9);
prb.solversettings = sdpsettings('solver','mosek','verbose',0,'mosek.MSK_DPAR_INTPNT_CO_TOL_PFEAS',1e-9,'mosek.MSK_DPAR_INTPNT_CO_TOL_REL_GAP',1e-9);

% prb.colors = {[1,0.8,0.5],[0.8,1,0.5],[0.5,1,0.8],[0.5,0.5,1]};
% prb.colors = {'magenta','red','blue','green'};
prb.colors = {[0,0.7,0.3],[1,0.65,0.3],[0.5,0.5,1],[1,0.5,0.5]};
prb.colors_trunk = {[0.1,0.1,0.2],[0.3,0.3,0.4],[0.5,0.5,0.6]};
prb.marker = {'o','s','^'};

prb.path_constraints = @(x,u) [norm(u) <= prb.umax;
                               u(3) >= prb.umin;
                               norm(u(1:2)) <= prb.tan_thet_tp * u(3)];

end