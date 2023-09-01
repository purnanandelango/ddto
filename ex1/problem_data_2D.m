function prb = problem_data_2D(K,scp_iters,wvc,wvb,wtr,cost_factor)

    prb.K = K;

    prb.ntarg = 4;    

    prb.n = 2;    

    prb.nx = 2*prb.n*prb.ntarg;
    prb.nu = (prb.n+1)*prb.ntarg;
    prb.np = 0;    

    % Convenient indices for state and input associated with trajectories to different targets
    [prb.idx_x,prb.idx_u,prb.idx_s] = util.disassemble_traj_idx(2*prb.n,prb.n,0,prb.ntarg);
    prb.idx_r = prb.idx_x(1:prb.n,:);
    prb.idx_v = prb.idx_x(prb.n+1:2*prb.n,:);    

    prb.tau = grid.generate_grid(0,1,K,'uniform');  % Generate grid in [0,1]

    prb.dtau = diff(prb.tau);
    
    prb.h = (1/10)*min(prb.dtau);                    % Step size for integration that computes FOH matrices
    prb.Kfine = 1+10*round(1/min(prb.dtau));         % Size of grid on which SCP solution is simulated

    % Deferrability index
    prb.Kstr = round(K/2);
    prb.taustr = prb.tau(prb.Kstr);
    [~,prb.Kstrfine] = min(abs(prb.taustr-grid.generate_grid(0,1,prb.Kfine,'uniform')));

    % System parameters

    prb.g = [0;0];                                  % External acceleration vector
        
    prb.c_d = 0.3;                                  % Drag coefficient
    
    % Bounds

    prb.rmax = 10;
    prb.vmax = 2.5;
    prb.umax = 4;
    prb.umin = 0.5;

    prb.smin     = 0.1;
    prb.smax     = 16;
    prb.dtmin    = 0.4;
    prb.dtmax    = 3;
    prb.ToFmax   = 50;
    prb.dumax    = 4;

    prb.snom = [1, 10];
    prb.ToFguess = 15;

    prb.cost_bound = 25*[1,1,1,1];

    % Obstacle avoidance
    prb.nobs = 2;

    prb.robs = [3  3;
                0  4];
    prb.aobs = [1.9  2];    
    
    % Boundary conditions
    prb.r1 = [0; 0];           
    prb.v1 = [0; 0];
    
    prb.rK = [2   5   0   5; 
              7   0   5   3];
              
    prb.vK = [1  0  0  0.5;
              0  1  0  1];

    assert(size(prb.rK,2) == prb.ntarg && size(prb.vK,2) == prb.ntarg,"Incorrect no. of targets states specified.");

    prb.x1 = reshape([prb.r1, prb.r1, prb.r1, prb.r1;
                      prb.v1, prb.v1, prb.v1, prb.v1],[prb.nx,1]);
    prb.xK = reshape([prb.rK;
                      prb.vK],[prb.nx,1]);
    prb.u1 = reshape([0.5*prb.umax*ones(prb.n,1), 0.5*prb.umax*ones(prb.n,1), 0.5*prb.umax*ones(prb.n,1), 0.5*prb.umax*ones(prb.n,1);
                      prb.ToFguess,               prb.ToFguess,               prb.ToFguess,               prb.ToFguess],[prb.nu,1]);
    prb.uK = reshape([0.5*prb.umax*ones(prb.n,1), 0.5*prb.umax*ones(prb.n,1), 0.5*prb.umax*ones(prb.n,1), 0.5*prb.umax*ones(prb.n,1);
                      prb.ToFguess,               prb.ToFguess,               prb.ToFguess,               prb.ToFguess],[prb.nu,1]);    

    % Scaling parameters

    % xmin = repmat([-0.5*prb.rmax*ones(prb.n,1); -0.5*prb.vmax*ones(prb.n,1)],[prb.ntarg,1]);
    % xmax = repmat([ 0.5*prb.rmax*ones(prb.n,1);  0.5*prb.vmax*ones(prb.n,1)],[prb.ntarg,1]);
    xmin = reshape([-0.5*prb.rmax*ones(prb.n,1),  -0.5*prb.rmax*ones(prb.n,1),  -0.5*prb.rmax*ones(prb.n,1),  -0.5*prb.rmax*ones(prb.n,1);
                    -0.5*prb.vmax*ones(prb.n,1),  -0.5*prb.vmax*ones(prb.n,1),  -0.5*prb.vmax*ones(prb.n,1),  -0.5*prb.vmax*ones(prb.n,1)],[prb.nx,1]);
    xmax = reshape([ 0.5*prb.rmax*ones(prb.n,1),   0.5*prb.rmax*ones(prb.n,1),   0.5*prb.rmax*ones(prb.n,1),   0.5*prb.rmax*ones(prb.n,1);
                     0.5*prb.vmax*ones(prb.n,1),   0.5*prb.vmax*ones(prb.n,1),   0.5*prb.vmax*ones(prb.n,1),   0.5*prb.vmax*ones(prb.n,1)],[prb.nx,1]);    
    
    % umin = repmat([prb.umin*ones(prb.n,1); prb.snom(1)],[prb.ntarg,1]);
    % umax = repmat([prb.umax*ones(prb.n,1); prb.snom(2)],[prb.ntarg,1]);
    umin = reshape([prb.umin*ones(prb.n,1), prb.umin*ones(prb.n,1), prb.umin*ones(prb.n,1), prb.umin*ones(prb.n,1);
                    prb.snom(1),            prb.snom(1),            prb.snom(1),            prb.snom(1)],[prb.nu,1]);
    umax = reshape([prb.umax*ones(prb.n,1), prb.umax*ones(prb.n,1), prb.umax*ones(prb.n,1), prb.umax*ones(prb.n,1);
                    prb.snom(2),            prb.snom(2),            prb.snom(2),            prb.snom(2)],[prb.nu,1]);

    [Sz,cz] = misc.generate_scaling({[xmin,xmax],[umin,umax]},[0,1]);

    prb.Sx = Sz{1}; prb.invSx = inv(Sz{1});
    prb.Su = Sz{2}; prb.invSu = inv(Sz{2});
    prb.cx = cz{1};
    prb.cu = cz{2};

    % SCP parameters

    prb.disc = "FOH";
    prb.foh_type = "v3_parallel";
    prb.scp_iters = scp_iters; % Maximum SCP iterations

    prb.solver_settings = sdpsettings('solver','gurobi','verbose',0,'gurobi.OptimalityTol',1e-9,'gurobi.FeasibilityTol',1e-9);
    
    prb.tr_norm = 2;
    
    prb.cost_term = @(x,u) norm(u);
    
    prb.wvc = wvc;
    prb.wvb = wvb; 
    prb.wtr = wtr;
    prb.cost_factor = cost_factor;
    
    prb.epsvc = 1e-8;
    prb.epstr = 1e-7;

    % Takes in unscaled data
    prb.time_of_maneuver = @(z,nu) disc.time_of_maneuver(prb.disc,prb.tau(1:prb.Kstr),nu(prb.n+1,1:prb.Kstr)); % Returns the time available to defer decision
    prb.time_grid = @(tau,z,s) disc.time_grid(prb.disc,tau,s);    
    
    % Convenient functions for accessing RHS of nonlinear and linearized ODE
    prb.dyn_func = @(t,z,nu)                    evaluate_dyn_func(z,nu,prb.ntarg,prb.n,prb.c_d,prb.g);
    prb.dyn_func_linearize = @(tbar,zbar,nubar) evaluate_linearization(zbar,nubar,prb.ntarg,prb.n,prb.c_d,prb.g);

end

function dz = evaluate_dyn_func(z,nu,ntarg,n,c_d,g)
    z = reshape(z,[2*n,ntarg]);
    nu = reshape(nu,[n+1,ntarg]);
    dz = zeros(2*n,ntarg);
    for j = 1:ntarg
        dz(:,j) = plant.doubleint.dyn_func(z(:,j),nu(1:n,j),nu(n+1,j),n,c_d,g);
    end
    dz = reshape(dz,[2*n*ntarg,1]);
end

function [A,B,w] = evaluate_linearization(z,nu,ntarg,n,c_d,g)
    z = reshape(z,[2*n,ntarg]);
    nu = reshape(nu,[n+1,ntarg]);
    A = [];
    B = [];
    w = zeros(2*n*ntarg,1);
    for j = 1:ntarg
        [Aj,Bj,Sj,wj] = plant.doubleint.compute_linearization(z(:,j),nu(1:n,j),nu(n+1,j),n,c_d,g);
        A = blkdiag(A,Aj);
        B = blkdiag(B,[Bj,Sj]);
        w((j-1)*(2*n)+1:j*2*n) = wj;
    end
end