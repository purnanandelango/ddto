function [cnstr,cost_fun,vc_cnstr] = sys_cnstr_cost(x,u,prb,...
                                                    ~,~)

    K = prb.K;
    n = prb.n;
    ntarg_p_1 = prb.ntarg+1;

    r = sdpvar(n,K,ntarg_p_1);
    v = sdpvar(n,K,ntarg_p_1);
    p = sdpvar(  K,ntarg_p_1);
    y = sdpvar(  K,ntarg_p_1);
    % T = sdpvar(n,K,ntarg_p_1);
    % s = sdpvar(  K,ntarg_p_1);   

    for k = 1:K
        xk = reshape(x(:,k),[2*n+2,ntarg_p_1]);
        % uk = reshape(u(:,k),[n+1,  ntarg_p_1]); 

        r(:,k,:) = xk(1:n,:);
        v(:,k,:) = xk(n+1:2*n,:);
        p(  k,:) = xk(2*n+1,:);
        y(  k,:) = xk(2*n+2,:);
        % T(:,k,:) = uk(1:n,:);
        % s(  k,:) = uk(n+1,:);
    end

    % Trunk initial conditions
    cnstr = [r(:,1,1) == prb.r1;
             v(:,1,1) == prb.v1;
             p(1,1) == prb.p1;
             y(1,1) == prb.y1];

    for k = 1:K
        cnstr = [cnstr;
                 prb.umin <= u(:,k) <= prb.umax];
    end

    for j=1:ntarg_p_1

        if j > 1
            % Branch final condition and trunk-branch stitching conditions            
            cnstr = [cnstr;
                     r(:,1,j) == r(:,K,1);
                     v(:,1,j) == v(:,K,1);
                     p(1,j) == p(K,1);
                     y(1,j) == y(K,1);
                     r(:,K,j) == prb.rK(:,j-1);
                     v(:,K,j) == prb.vK(:,j-1);
                     ];

            % Cumulative trajectory cost constraint            
            cnstr = [cnstr;
                     p(K,j) <= prb.cost_bound(j-1)];
        
        end

        for k = 1:K-1              
            cnstr = [cnstr;
                     y(k+1,j) - y(k,j) <= prb.eps_cnstr];           
        end

    end

    % Scaled dilation factors of trunk trajectory
    s1_scl = prb.invSu(n+1,n+1)*(u(n+1,:) - prb.cu(n+1));

    % Maximize duration of trunk trajectory
    cost_fun = -prb.cost_factor*sum(s1_scl(1:end-1));

    vc_cnstr = 0;

end