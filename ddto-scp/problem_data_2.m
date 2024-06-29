function prb = problem_data_2(K,scp_iters,w_ep,w_px,cost_factor,...
                            zinit,rvtarg,cost_bound)
    
    prb.K = K;

    % Dimension of system

    prb.ntarg = size(rvtarg,2); ntarg = prb.ntarg; % No. of targets

    % prb.m = 7; % No. of path constraints

    prb.n = 3; n = prb.n;

    prb.nx = (2*n+1+1)*(ntarg+1);
    prb.nu = (n+1)*(ntarg+1);
    
    % Generate grid in [0,1]
    prb.tau = grid.generate_grid(0,1,K,'uniform'); 
    prb.dtau = diff(prb.tau); min_dtau = min(prb.dtau);
    
    prb.h = (1/19)*min_dtau;                    % Step size for integration that computes discretization matrices
    prb.Kfine = 1+100*round(1/min_dtau);         % Size of grid on which SCP solution is simulated
    
    % System parameters

    prb.accl        = [0;0;-9.806];                 % External acceleration    
    prb.c_d         = 0.01;                     % Drag coefficient
    
    % Bounds

    prb.rmax        = 40;
    prb.vmax        = 08;
    prb.pmin        = 0;
    prb.pmax        = 1000;

    prb.cost_bound  = cost_bound;

    prb.ymin        = 0;
    prb.ymax        = 1;

    prb.Tmin        = 5.0;
    prb.Tmax        = 20.0;    

    prb.ehat        = [0;0;1];
    prb.deltamax    = 60;                       % [deg] 

    prb.smin        = 1;
    prb.smax        = 15;
    prb.ToFguess    = 10;

    % Obstacle avoidance

    prb.nobs        = 2;

    % Centers
    prb.qobs        = [-5 -10;                  
                        1  20;
                       10  10];

    % Shape matrices
    prb.robs        = {[0.2,0.1,0.2],[0.1,0.2,0.2]};

    prb.Hobs        = {blkdiag(diag(prb.robs{1}(1:2))*geom.rot_mat_2D(0),prb.robs{1}(3)), ...
                       blkdiag(diag(prb.robs{2}(1:2))*geom.rot_mat_2D(0),prb.robs{1}(3))}; 

 
    % Boundary conditions

    prb.r1 = zinit(1:n);
    prb.v1 = zinit(n+1:2*n);
    prb.p1 = zinit(2*n+1);
    prb.y1 = zinit(2*n+2);

    prb.rK = rvtarg(1:n,:);
    prb.vK = rvtarg(n+1:2*n,:);

    prb.rmid = sum(prb.rK + repmat(prb.r1,[1,ntarg]),2)/(2*ntarg);
    prb.vmid = sum(prb.vK + repmat(prb.v1,[1,ntarg]),2)/(2*ntarg);
    % prb.Tmid = 0.5*(prb.Tmax+prb.Tmin)*ones(n,1);
    prb.Tmid = prb.Tmin*ones(n,1);

    % Initialization generator

    prb.x1 = reshape([[prb.r1;
                       prb.v1;
                       prb.p1;
                       prb.y1], repmat([prb.rmid;prb.vmid;prb.p1;prb.y1],[1,ntarg])],[prb.nx,1]);

    prb.xK = reshape([[prb.rmid, prb.rK;
                       prb.vmid, prb.vK];
                      repmat([prb.pmax;prb.y1],[1,ntarg+1])],[prb.nx,1]);    

    prb.u1 = reshape(repmat([prb.Tmid;
                             prb.ToFguess],[1,ntarg+1]),[prb.nu,1]);
    prb.uK = reshape(repmat([prb.Tmid;
                             prb.ToFguess],[1,ntarg+1]),[prb.nu,1]);

    % Scaling parameters
    xmin = reshape(repmat([-0.5*prb.rmax*ones(prb.n,1);
                           -0.5*prb.vmax*ones(prb.n,1);
                                prb.pmin; 
                                prb.ymin],[1,ntarg+1]),[prb.nx,1]);
    xmax = reshape(repmat([0.5*prb.rmax*ones(prb.n,1);
                           0.5*prb.vmax*ones(prb.n,1);
                               prb.pmax; 
                               prb.ymax],[1,ntarg+1]),[prb.nx,1]);    
    
    umin = reshape(repmat([-prb.Tmax*ones(prb.n,1);
                            prb.smin],[1,ntarg+1]),[prb.nu,1]);     prb.umin = umin;
    umax = reshape(repmat([ prb.Tmax*ones(prb.n,1);
                            prb.smax],[1,ntarg+1]),[prb.nu,1]);     prb.umax = umax;

    [Sz,cz] = misc.generate_scaling({[xmin,xmax],[umin,umax]},[0,1]);

    prb.Sx = Sz{1}; prb.invSx = inv(Sz{1});
    prb.Su = Sz{2}; prb.invSu = inv(Sz{2});
    prb.cx = cz{1};
    prb.cu = cz{2};

    prb.eps_cnstr = 1e-5;

    cnstr_scl = diag([...
                      1;
                      1;
                      0.1;
                      0.1; 
                      0.1;
                      0.1;
                      1;
                      ]);
    cnstr_buffer = [0;
                    0;
                    0;
                    0;
                    0;
                    0;
                    0;
                    ];

    % Constraint parameters

    prb.cnstr_fun       = @(rvp,T) cnstr_scl*[-norm(prb.Hobs{1}*(rvp(1:n)-prb.qobs(:,1)))^2 + 1;            % Obstacle avoidance
                                              -norm(prb.Hobs{2}*(rvp(1:n)-prb.qobs(:,2)))^2 + 1;            % Obstacle avoidance
                                               norm(rvp(n+1:2*n))^2 - prb.vmax^2;                           % Speed upperbound
                                               norm(T)^2 - prb.Tmax^2;                                      % Thrust upperbound
                                              -norm(T)^2 + prb.Tmin^2;                                      % Thrust lowerbound
                                               norm(T)^2 - (secd(prb.deltamax)*prb.ehat'*T)^2;              % Thrust pointing
                                               -prb.ehat'*T
                                             ] ... 
                                             + cnstr_buffer;

    prb.cnstr_fun_jac_rvp = @(rvp,T) cnstr_scl*[-2*(rvp(1:n)-prb.qobs(:,1))'*prb.Hobs{1}'*prb.Hobs{1}, zeros(1,n),      0;
                                                -2*(rvp(1:n)-prb.qobs(:,2))'*prb.Hobs{2}'*prb.Hobs{2}, zeros(1,n),      0;
                                                 zeros(1,n),                                           2*rvp(n+1:2*n)', 0;
                                                 zeros(1,n),                                           zeros(1,n),      zeros(1,1);
                                                 zeros(1,n),                                           zeros(1,n),      zeros(1,1);
                                                 zeros(1,n),                                           zeros(1,n),      zeros(1,1);
                                                 zeros(1,n),                                           zeros(1,n),      zeros(1,1);
                                             ];

    prb.cnstr_fun_jac_T = @(rvp,T) cnstr_scl*[ zeros(3,n);
                                               2*T';
                                              -2*T';
                                               2*T' - 2*(secd(prb.deltamax)^2)*(prb.ehat'*T)*prb.ehat';
                                               -prb.ehat';
                                             ];

    % SCP parameters

    prb.disc = "ZOH";
    prb.zoh_type = "v3_parallel";

    prb.ode_solver = {'ode45',odeset('RelTol',1e-5,'AbsTol',1e-7)};
    prb.scp_iters = scp_iters; % Maximum SCP iterations

    % prb.solver_settings = sdpsettings('solver','quadprog','verbose',false);    
    % prb.solver_settings = sdpsettings('solver','ecos','verbose',false,'ecos.abstol',1e-8,'ecos.reltol',1e-8);    
    prb.solver_settings = sdpsettings('solver','gurobi','verbose',false,'gurobi.OptimalityTol',1e-9,'gurobi.FeasibilityTol',1e-9);
    % prb.solver_settings = sdpsettings('solver','mosek','verbose',false,'mosek.MSK_DPAR_INTPNT_CO_TOL_PFEAS',1e-9,'mosek.MSK_DPAR_INTPNT_CO_TOL_REL_GAP',1e-9);
    % prb.solver_settings = sdpsettings('solver','osqp','verbose',false,'osqp.eps_abs',1e-9,'osqp.eps_rel',1e-9,'osqp.max_iter',1e5);        
    
    % prb.px_norm = 2;
    % prb.px_norm = inf;   
    prb.px_norm = 'quad';
    
    prb.w_ep = w_ep;
    prb.w_px = w_px;
    prb.cost_factor = cost_factor;
    
    prb.eps_ep = 1e-7;
    prb.eps_px = 1e-3;

    % Time of maneuver and time grid
    prb.time_of_maneuver =     @(x,u) disc.time_of_maneuver(prb.disc,prb.tau,u(n+1,:)); % Returns the time available to defer decision
    prb.time_grid        = @(tau,x,u)        disc.time_grid(prb.disc,    tau,u(n+1,:));    
    
    % Convenient functions for accessing RHS of nonlinear and linearized ODE
    prb.dyn_func           = @(tau,x,u)   evaluate_dyn_func     (x,u,n,ntarg+1,prb.c_d,prb.accl,prb.cnstr_fun);
    prb.dyn_func_linearize = @(tau,x,u)   evaluate_linearization(x,u,n,ntarg+1,prb.c_d,prb.accl,prb.cnstr_fun,...
                                                                                       prb.cnstr_fun_jac_rvp,prb.cnstr_fun_jac_T);

end

function f = evaluate_dyn_func(x,u,n,ntarg_p_1,c_d,accl,cnstr_fun)

    x = reshape(x,[2*n+2,ntarg_p_1]);
    u = reshape(u,[n+1,ntarg_p_1]);

    rvp = x(1:2*n+1,1:ntarg_p_1);
    v = x(n+1:2*n,1:ntarg_p_1);
    T = u(1:n,1:ntarg_p_1);
    s = u(n+1,1:ntarg_p_1);

    f = zeros(2*n+2,ntarg_p_1);
    for j = 1:ntarg_p_1
        cnstr_val = cnstr_fun(rvp(:,j),T(:,j));
        
        F = [v(:,j);
             T(:,j) + accl - c_d*norm(v(:,j))*v(:,j);
             norm(T(:,j))^2];
    
        f(:,j) = s(j)*[F;
                       sum( arrayfun(@(y) max(0,y)^2, cnstr_val) )];
    end

    f = f(:);
end

function [A,B,w] = evaluate_linearization(x,u,n,ntarg_p_1,c_d,accl,cnstr_fun, ...
                                                              cnstr_fun_jac_rvp,cnstr_fun_jac_T)

    x = reshape(x,[2*n+2,ntarg_p_1]);
    u = reshape(u,[n+1,ntarg_p_1]);

    rvp = x(1:2*n+1,1:ntarg_p_1);
    v = x(n+1:2*n,1:ntarg_p_1);
    T = u(1:n,1:ntarg_p_1);
    s = u(n+1,1:ntarg_p_1);

    A = [];
    B = [];

    for j = 1:ntarg_p_1

        cnstr_val         = cnstr_fun(rvp(:,j),T(:,j));
        abs_cnstr_val     = arrayfun(@(y) max(0,y),cnstr_val);
        abs2_cnstr_val    = arrayfun(@(y) max(0,y)^2,cnstr_val);  
        cnstr_val_jac_rvp = cnstr_fun_jac_rvp(rvp(:,j),T(:,j));
        cnstr_val_jac_T   = cnstr_fun_jac_T(rvp(:,j),T(:,j));
    
        if norm(v(:,j)) < 1e-7
            term_v = 0;
        else
            term_v = v(:,j)*v(:,j)'/norm(v(:,j)); 
        end
    
        F = [v(:,j);
             T(:,j) + accl - c_d*norm(v(:,j))*v(:,j);
             norm(T(:,j))^2];
    
        dFdrvp = [zeros(n),  eye(n),                           zeros(n,1);
                  zeros(n), -c_d*(norm(v(:,j))*eye(n)+term_v), zeros(n,1);
                  zeros(1,2*n+1)];
    
        dFdT = [zeros(n);
                eye(n);
                2*T(:,j)'];
    
        A = blkdiag(A,s(j)*[dFdrvp,                             zeros(2*n+1,1);
                            2*abs_cnstr_val'*cnstr_val_jac_rvp, 0]);
        
        B = blkdiag(B,[s(j)*dFdT,                              F;
                       2*s(j)*abs_cnstr_val'*cnstr_val_jac_T,  sum(abs2_cnstr_val)]);        

    end
    
    f = evaluate_dyn_func(x(:),u(:),n,ntarg_p_1,c_d,accl,cnstr_fun);
    
    w = f - A*x(:) - B*u(:);
end