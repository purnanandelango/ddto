function [cnstr,cost_fun,vcvb_cnstr] = sys_cnstr_cost(z,nu,prb,...
                                                      zbar,nubar)
    K = prb.K;
    n = prb.n;
    ntarg = prb.ntarg;

    % Containers for unscaled variables
    r   = sdpvar(n,K,ntarg);
    v   = sdpvar(n,K,ntarg);
    u   = sdpvar(n,K,ntarg);
    s   = sdpvar(K,ntarg);
    dt  = sdpvar(K-1,ntarg); % Time elapsed in each interval

    % Obstacle avoidance buffer
    nu_ncvx = sdpvar(prb.nobs+1,ntarg,K);    

    rbar = zeros(n,K,ntarg);
    vbar = zeros(n,K,ntarg);
    ubar = zeros(n,K,ntarg);

    % Define unscaled states and control inputs
    for j = 1:ntarg
        for k = 1:K
            r(:,k,j)   = prb.Sx(prb.idx_r(:,j),prb.idx_r(:,j))     *z(prb.idx_r(:,j),k)    + prb.cx(prb.idx_r(:,j));
            v(:,k,j)   = prb.Sx(prb.idx_v(:,j),prb.idx_v(:,j))     *z(prb.idx_v(:,j),k)    + prb.cx(prb.idx_v(:,j));
    
            u(:,k,j)   = prb.Su(prb.idx_u(:,j),prb.idx_u(:,j))     *nu(prb.idx_u(:,j),k)   + prb.cu(prb.idx_u(:,j));        
            s(k,j)     = prb.Su(prb.idx_s(j),prb.idx_s(j))         *nu(prb.idx_s(j),k)     + prb.cu(prb.idx_s(j));

            rbar(:,k,j) = zbar(prb.idx_r(:,j),k);
            vbar(:,k,j) = zbar(prb.idx_v(:,j),k);  
            ubar(:,k,j) = nubar(prb.idx_u(:,j),k);
        end
    end    

    cnstr = [];

    for j = 1:ntarg

        % Boundary conditions
        cnstr = [cnstr;
                 r(:,1,j)   == prb.r1;
                 v(:,1,j)   == prb.v1;
                 r(:,K,j)   == prb.rK(:,j);
                 v(:,K,j)   == prb.vK(:,j)];       

        % Constraints
        for k = 1:K   

            if k < K
                switch prb.disc
                    case "ZOH"
                        dt(k,j) = prb.dtau(k)*s(k,j); 
                    case "FOH"
                        dt(k,j) = 0.5*prb.dtau(k)*(s(k+1,j)+s(k,j));
                end
                % Upper and lower bound on time elapsed in each interval
                cnstr = [cnstr; 
                         prb.dtmin <= dt(k,j) <= prb.dtmax];
                cnstr = [cnstr;      
                         norm(u(:,k+1,j)-u(:,k,j)) <= dt(k,j)*prb.dumax];   % Upper bound on magnitude of rate of change of thrust
            end            
            
            cnstr = [cnstr;
                     norm(u(:,k,j)) <= prb.umax;                            % Thrust upper bound                                                                   % Thrust magnitude upper bound
                     norm(v(:,k,j)) <= prb.vmax;                            % Velocity magnitude upper bound
                     -prb.rmax <= r(:,k,j) <= prb.rmax;                     % Bounds on position 
                     prb.smin <= s(k,j) <= prb.smax];                       % Lower and upper bounds on dilation factor                                                         % Lower and upper bounds on dilation factor

            % Deferrability
            if j > 1 && k <= prb.Kstr
                cnstr = [cnstr;
                         r(:,k,1) == r(:,k,j);
                         v(:,k,1) == v(:,k,j)];
            end

            % stage_cost(k,j) = prb.cost_term(z(prb.idx_x(:,j),k),nu(prb.idx_u(:,j),k));      % Scaled
            stage_cost(k,j) = prb.cost_term([r(:,k,j);v(:,k,j)],u(:,k,j));                  % Unscaled                        

            for i = 1:prb.nobs
                cnstr = [cnstr;
                         norm(rbar(:,k,j)-prb.robs(:,i)) - prb.aobs(i) + dot(rbar(:,k,j)-prb.robs(:,i),r(:,k,j)-rbar(:,k,j))/norm(rbar(:,k,j)-prb.robs(:,i)) + nu_ncvx(i,j,k) >= 0;
                         nu_ncvx(i,j,k) >= 0];
            end  

            % Thrust lower bound
            cnstr = [cnstr;
                     norm(ubar(:,k,j)) - prb.umin + dot(ubar(:,k,j),u(:,k,j)-ubar(:,k,j))/norm(ubar(:,k,j)) + nu_ncvx(end,j,k) >= 0;
                     nu_ncvx(end,j,k) >= 0];
        
        end
    
        % Suboptimality constraint
        cnstr = [cnstr; sum(stage_cost(:,j)) <= prb.cost_bound(j)];    

        % Time of maneuver upper bound
        cnstr = [cnstr; sum(dt(:,j)) <= prb.ToFmax];                              

    end  

    % Time available to defer decision
    defer_time = sum(dt(1:prb.Kstr,1));

    vcvb_cnstr = sum(nu_ncvx(:));    

    % Maximize deferrability
    cost_fun = -prb.cost_factor*defer_time + prb.wvb*vcvb_cnstr;        

end