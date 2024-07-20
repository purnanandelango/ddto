function sol = solve_ddto_micp_relax_v1(prb,weights)

    nx  = prb.nx;
    nu  = prb.nu;
    n   = prb.n;
    N   = prb.N;
    i   = prb.i;
    Ni  = prb.Ni;

    yalmip clear

    X   = cell(1,n);
    U   = cell(1,n);
    for j=1:n
        X{j} = sdpvar(nx,N(j));
        U{j} = sdpvar(nu,N(j)-1);
    end

    cnstr = [];
    
    for j = 1:n
        cum_cost = 0;
        for k = 1:N(j)-1
            % Dynamics constraint, state & input constraints
            cnstr = [cnstr;
                     X{j}(:,k+1) == prb.A*X{j}(:,k) + prb.B*U{j}(:,k) + prb.c;
                     norm(U{j}(:,k)) <= prb.umax;
                     U{j}(3,k) >= prb.umin;
                     norm(U{j}(1:2,k)) <= prb.tan_thet_tp*U{j}(3,k);
                     ];
            cum_cost = cum_cost + prb.stage_cost(X{j}(:,k),U{j}(:,k));
        end
        cnstr = [cnstr;
                 cum_cost <= prb.lmax;
                 X{j}(:,1) == prb.z0;
                 X{j}(:,N(j)) == prb.zf(:,j);                 
                 ];
    end

    objfun = 0;

    if isempty(weights)
        weights = cell(1,n);
        for j = setdiff(1:n,i)
            weights{j} = ones(Ni(j),1);
        end        
    end

    Xi  = cell(1,n); 
    for j = setdiff(1:n,i)
        Xi{j} = sdpvar(Ni(j),1);
    end    

    for j = setdiff(1:n,i)    
        for k = 1:Ni(j)
            cnstr = [cnstr;
                     norm(X{j}(:,k) - X{i}(:,k)) <= Xi{j}(k)*prb.M;
                     0 <= Xi{j}(k) <= 1;
                     ];
            objfun = objfun + Xi{j}(k) * weights{j}(k);
        end
    end

    yalmip_out = optimize(cnstr,objfun,prb.solversettings);

    sol = struct;
    sol.X   = cell(1,n);
    sol.U   = cell(1,n);    
    sol.Xi  = cell(1,n);    
    for j=1:n
        sol.X{j} = value(X{j});
        sol.U{j} = value(U{j});
        if j ~= i
            sol.Xi{j} = value(Xi{j});        
        end
    end

end