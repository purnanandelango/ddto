function prb = problem_data_2D(K,scp_iters,wvc,wtr,cost_factor)

    prb.K = K;

    prb.ntarg = 4;    

    prb.n = 2;    

    prb.m = 4;

    prb.nx = (2*prb.n + prb.m)*prb.ntarg;
    prb.nu = (prb.n + 1)*prb.ntarg;
    prb.np = 0;    

    % Convenient indices for state and input associated with trajectories to different targets
    prb.idx_r   = zeros(prb.n,prb.ntarg);
    prb.idx_v   = zeros(prb.n,prb.ntarg);
    prb.idx_bet = zeros(prb.m,prb.ntarg);
    prb.idx_T   = zeros(prb.n,prb.ntarg);
    prb.idx_s   = zeros(1,prb.ntarg);
    for j = 1:prb.ntarg
        prb.idx_r(:,j)   = (j-1)*(2*prb.n+prb.m)+1         : (j-1)*(2*prb.n+prb.m)+prb.n;
        prb.idx_v(:,j)   = (j-1)*(2*prb.n+prb.m)+prb.n+1   : (j-1)*(2*prb.n+prb.m)+2*prb.n;
        prb.idx_bet(:,j) = (j-1)*(2*prb.n+prb.m)+2*prb.n+1 :  j*(2*prb.n+prb.m);
        prb.idx_T(:,j)   = (j-1)*(prb.n+1)+1               : (j-1)*(prb.n+1)+prb.n;
        prb.idx_s(j)     = j*(prb.n+1);    
    end    

    prb.tau = grid.generate_grid(0,1,K,'uniform');  % Generate grid in [0,1]

    prb.dtau = diff(prb.tau);
    
    prb.h = (1/40)*min(prb.dtau);                    % Step size for integration that computes FOH matrices
    prb.Kfine = 1+round(20/min(prb.dtau));           % Size of grid on which SCP solution is simulated

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
    prb.dTmax    = 4;

    prb.betmin = 0.01*ones(prb.m,1);
    prb.betmax = 0.1*[1;
                  1;
                  1;
                  1];    

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
    prb.bet1 = zeros(prb.m,1);
    
    prb.rK = [2   5   0   5; 
              7   0   5   3];
              
    prb.vK = [1  0  0  0.5;
              0  1  0  1];

    prb.betK = zeros(prb.m,1);    

    assert(size(prb.rK,2) == prb.ntarg && size(prb.vK,2) == prb.ntarg,"Incorrect no. of targets states specified.");

    prb.x1 = reshape([prb.r1,   prb.r1,   prb.r1,   prb.r1;
                      prb.v1,   prb.v1,   prb.v1,   prb.v1;
                      prb.bet1, prb.bet1, prb.bet1, prb.bet1],[prb.nx,1]);
    prb.xK = reshape([prb.rK;
                      prb.vK; ...
                      prb.betK, prb.betK, prb.betK, prb.betK],[prb.nx,1]);
    prb.u1 = reshape([0.5*prb.umax*ones(prb.n,1), 0.5*prb.umax*ones(prb.n,1), 0.5*prb.umax*ones(prb.n,1), 0.5*prb.umax*ones(prb.n,1);
                      prb.ToFguess,               prb.ToFguess,               prb.ToFguess,               prb.ToFguess],[prb.nu,1]);
    prb.uK = reshape([0.5*prb.umax*ones(prb.n,1), 0.5*prb.umax*ones(prb.n,1), 0.5*prb.umax*ones(prb.n,1), 0.5*prb.umax*ones(prb.n,1);
                      prb.ToFguess,               prb.ToFguess,               prb.ToFguess,               prb.ToFguess],[prb.nu,1]);    

    % Scaling parameters

    % xmin = repmat([-0.5*prb.rmax*ones(prb.n,1); -0.5*prb.vmax*ones(prb.n,1); prb.betmin],[prb.ntarg,1]);
    % xmax = repmat([ 0.5*prb.rmax*ones(prb.n,1);  0.5*prb.vmax*ones(prb.n,1); prb.betmax],[prb.ntarg,1]);
    xmin = reshape([-0.5*prb.rmax*ones(prb.n,1),  -0.5*prb.rmax*ones(prb.n,1),  -0.5*prb.rmax*ones(prb.n,1),  -0.5*prb.rmax*ones(prb.n,1);
                    -0.5*prb.vmax*ones(prb.n,1),  -0.5*prb.vmax*ones(prb.n,1),  -0.5*prb.vmax*ones(prb.n,1),  -0.5*prb.vmax*ones(prb.n,1);
                     prb.betmin,                   prb.betmin,                   prb.betmin,                   prb.betmin],[prb.nx,1]);
    xmax = reshape([ 0.5*prb.rmax*ones(prb.n,1),   0.5*prb.rmax*ones(prb.n,1),   0.5*prb.rmax*ones(prb.n,1),   0.5*prb.rmax*ones(prb.n,1);
                     0.5*prb.vmax*ones(prb.n,1),   0.5*prb.vmax*ones(prb.n,1),   0.5*prb.vmax*ones(prb.n,1),   0.5*prb.vmax*ones(prb.n,1);
                     prb.betmax,                   prb.betmax,                   prb.betmax,                   prb.betmax],[prb.nx,1]);    
    
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
    prb.foh_type = "v3";
    prb.scp_iters = scp_iters; % Maximum SCP iterations

    prb.solver_settings = sdpsettings('solver','gurobi','verbose',0);
    
    prb.tr_norm = 2;
    
    prb.cost_term = @(z) norm(z);
    prb.subopt_type = 'sum_stage_cost';
    
    % prb.subopt_type = 'sum_quad_u';
    
    prb.wvc = wvc;
    % prb.wvb = wvb; 
    prb.wtr = wtr;
    prb.cost_factor = cost_factor;
    
    prb.epsvc = 1e-8;
    prb.epstr = 1e-7;

    % Takes in unscaled data
    prb.time_of_maneuver = @(x,u) disc.time_of_maneuver(prb.disc,prb.tau(1:prb.Kstr),u(prb.n+1,1:prb.Kstr)); % Returns the time available to defer decision
    prb.time_grid = @(tau,x,s) disc.time_grid(prb.disc,tau,s);    
    
    % Convenient functions for accessing RHS of nonlinear and linearized ODE
    prb.dyn_func = @(t,x,u)                    evaluate_dyn_func(x,u,           prb.ntarg,prb.n,prb.m,prb.c_d,prb.g,prb.robs,prb.aobs,prb.vmax,prb.umin);
    prb.dyn_func_linearize = @(tbar,xbar,ubar) evaluate_linearization(xbar,ubar,prb.ntarg,prb.n,prb.m,prb.c_d,prb.g,prb.robs,prb.aobs,prb.vmax,prb.umin);

end

function dx = evaluate_dyn_func(x,u,ntarg,n,m,c_d,g,robs,aobs,vmax,umin)
    x = reshape(x,[2*n+m,ntarg]);
    u = reshape(u,[n+1,ntarg]);
    dx = zeros(2*n+m,ntarg);
    for j = 1:ntarg
        dx(:,j) = [plant.doubleint.dyn_func(x(1:2*n,j),u(1:n,j),u(n+1,j),n,c_d,g);
                   max(0, -norm(x(1:n,j)-robs(:,1)) + aobs(1) )^2;
                   max(0, -norm(x(1:n,j)-robs(:,2)) + aobs(2) )^2;
                   max(0,  norm(x(n+1:2*n,j))^2 - vmax^2      )^2;
                   max(0, -norm(u(1:n,j)) + umin              )^2]; 
    end
    dx = reshape(dx,[(2*n+m)*ntarg,1]);
end

function [A,B,w] = evaluate_linearization(xvec,uvec,ntarg,n,m,c_d,g,robs,aobs,vmax,umin)
    x = reshape(xvec,[2*n+m,ntarg]);
    u = reshape(uvec,[n+1,ntarg]);
    f = evaluate_dyn_func(x,u,ntarg,n,m,c_d,g,robs,aobs,vmax,umin);
    
    A = [];
    B = [];
    for j = 1:ntarg
        [Aj,Bj,Sj] = plant.doubleint.compute_linearization(x(1:2*n,j),u(1:n,j),u(n+1,j),n,c_d,g);
        
        cnstr = [-norm(x(1:n,j)-robs(:,1)) + aobs(1);
                 -norm(x(1:n,j)-robs(:,2)) + aobs(2);
                  norm(x(n+1:2*n,j))^2 - vmax^2;
                 -norm(u(1:n,j)) + umin];
        cnstr_jac_state = [-(x(1:n,j)-robs(:,1))'/norm(x(1:n,j)-robs(:,1)) zeros(1,n);
                           -(x(1:n,j)-robs(:,2))'/norm(x(1:n,j)-robs(:,2)) zeros(1,n);
                            zeros(1,n) 2*x(n+1:2*n,j)'
                            zeros(1,2*n)];
        cnstr_jac_ctrl = [ zeros(1,n+1);
                           zeros(1,n+1);
                           zeros(1,n+1);
                          -u(1:n,j)'/norm(u(1:n,j)) 0];
        
        A_cnstr = [2*max(cnstr(1),0)*cnstr_jac_state(1,:);
                   2*max(cnstr(2),0)*cnstr_jac_state(2,:);
                   2*max(cnstr(3),0)*cnstr_jac_state(3,:);
                   2*max(cnstr(4),0)*cnstr_jac_state(4,:)];

        B_cnstr = [2*max(cnstr(1),0)*cnstr_jac_ctrl(1,:);
                   2*max(cnstr(2),0)*cnstr_jac_ctrl(2,:);
                   2*max(cnstr(3),0)*cnstr_jac_ctrl(3,:);
                   2*max(cnstr(4),0)*cnstr_jac_ctrl(4,:)];
        
        A = blkdiag(A,[Aj,      zeros(2*n,m);
                       A_cnstr, zeros(m,m)]);

        B = blkdiag(B,[Bj,Sj;
                       B_cnstr]);
                       
    end
    w = f - A*xvec - B*uvec;
end