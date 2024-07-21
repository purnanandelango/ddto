function sol = solve_ddto_micp(prb)

    nx  = prb.nx;
    nu  = prb.nu;
    n   = prb.n;
    N   = prb.N;
    i   = prb.i;
    Ni  = prb.Ni;

    yalmip clear

    X   = cell(1,n);
    U   = cell(1,n);
    for j = 1:n
        X{j} = sdpvar(nx,N(j));
        U{j} = sdpvar(nu,N(j)-1);
    end

    cnstr = [];
    
    cum_cost_var = sdpvar(1,n);
    for j = 1:n
        cum_cost = 0;
        for k = 1:N(j)-1
            % Dynamics constraint and path constraints
            cnstr = [cnstr;
                     X{j}(:,k+1) == prb.A*X{j}(:,k) + prb.B*U{j}(:,k) + prb.c;
                     prb.path_constraints(X{j}(:,k),U{j}(:,k));
                     ];
            % Cumulative trajectory cost constraint
            cum_cost = cum_cost + prb.stage_cost(X{j}(:,k),U{j}(:,k));
        end
        cum_cost_var(j) = cum_cost;        
        cnstr = [cnstr;
                 cum_cost_var(j) <= prb.lmax;
                 X{j}(:,1) == prb.z0;
                 X{j}(:,N(j)) == prb.zf(:,j);                 
                 ];
    end

    % objfun = 0.01*sum(cum_cost_var);
    objfun = 0;

    Xi  = cell(1,n);    
    for j = setdiff(1:n,i)
        Xi{j} = binvar(Ni(j),1);
    end    

    for j = setdiff(1:n,i)    
        for k = 1:Ni(j)
            cnstr = [cnstr;
                     norm(X{j}(:,k) - X{i}(:,k)) <= Xi{j}(k)*prb.M;
                     ];
            objfun = objfun + Xi{j}(k);
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
    % sol.objval = value(objfun);
    sol.objval = compute_ddto_cost(sol.X,i);
    sol.cum_cost = value(cum_cost_var); 
    sol.parsetime = yalmip_out.yalmiptime;
    sol.solvetime = yalmip_out.solvertime;

end