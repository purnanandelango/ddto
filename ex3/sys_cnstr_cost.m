function [cnstr,cost_fun,vcvb_cnstr] = sys_cnstr_cost(z,nu,prb,...
                                                      zbar,nubar)
    K = prb.K;
    n = prb.n;
    m = prb.m;
    ntarg = prb.ntarg;

    % Containers for unscaled variables
    r   = sdpvar(n,K,ntarg+1);
    v   = sdpvar(n,K,ntarg+1);
    bet = sdpvar(m,K,ntarg+1);
    u   = sdpvar(n,K,ntarg+1);
    s   = sdpvar(K,ntarg+1);
    dt  = sdpvar(K-1,ntarg+1); % Time elapsed in each interval    

    stage_cost = sdpvar(K,ntarg+1);    

    rbar = zeros(n,K,ntarg+1);
    vbar = zeros(n,K,ntarg+1);
    ubar = zeros(n,K,ntarg+1);

    % Define unscaled states and control inputs
    for j = 1:ntarg+1
        for k = 1:K
            r(:,k,j)   = prb.Sx(prb.idx_r(:,j),prb.idx_r(:,j))     *z(prb.idx_r(:,j),k)    + prb.cx(prb.idx_r(:,j));
            v(:,k,j)   = prb.Sx(prb.idx_v(:,j),prb.idx_v(:,j))     *z(prb.idx_v(:,j),k)    + prb.cx(prb.idx_v(:,j));
            bet(:,k,j) = prb.Sx(prb.idx_bet(:,j),prb.idx_bet(:,j)) *z(prb.idx_bet(:,j),k)  + prb.cx(prb.idx_bet(:,j));
    
            u(:,k,j)   = prb.Su(prb.idx_u(:,j),prb.idx_u(:,j))     *nu(prb.idx_u(:,j),k)   + prb.cu(prb.idx_u(:,j));        
            s(k,j)     = prb.Su(prb.idx_s(j),prb.idx_s(j))         *nu(prb.idx_s(j),k)     + prb.cu(prb.idx_s(j));

            rbar(:,k,j) = zbar(prb.idx_r(:,j),k);
            vbar(:,k,j) = zbar(prb.idx_v(:,j),k);  
            ubar(:,k,j) = nubar(prb.idx_u(:,j),k);
        end
    end    

    cnstr = [];    

    for j = 1:ntarg+1

        % Boundary conditions
        if j == 1
            cnstr = [cnstr;
                     r(:,1,1)   == prb.r1;
                     v(:,1,1)   == prb.v1];            
        else
            cnstr = [cnstr;
                     r(:,1,j)   == r(:,K,1);
                     v(:,1,j)   == r(:,K,1);
                     u(:,1,j)   == u(:,K,1);
                     s(1,j)     == s(K,1);
                     r(:,K,j)   == prb.rK(:,j-1);
                     v(:,K,j)   == prb.vK(:,j-1)];       
        end

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
            
            % Convex constraints on control
            cnstr = [cnstr;
                     norm(u(:,k,j)) <= prb.umax;        % Thrust upper bound
                     prb.smin <= s(k,j) <= prb.smax];   % Lower and upper bounds on dilation factor

            % Integrated constraint violation
            cnstr = [cnstr;
                     bet(:,k,j) == zeros(m,1)];

            % stage_cost(k,j) = prb.cost_term(z(prb.idx_x(:,j),k),nu(prb.idx_u(:,j),k));      % Scaled
            stage_cost(k,j) = prb.cost_term([r(:,k,j);v(:,k,j)],u(:,k,j));                  % Unscaled                        
        
        end

        if j > 1
            % Suboptimality constraint            
            cnstr = [cnstr; sum(stage_cost(:,1)) + sum(stage_cost(:,j)) <= prb.cost_bound(j-1)];            
            % Time of maneuver upper bound
            cnstr = [cnstr; sum(dt(:,1)) + sum(dt(:,j)) <= prb.ToFmax];            
        end                                        

    end

    % Time available to defer decision
    defer_time = sum(dt(:,1));

    vcvb_cnstr = 0;

    % Maximize deferrability
    cost_fun = -prb.cost_factor*defer_time;        

end