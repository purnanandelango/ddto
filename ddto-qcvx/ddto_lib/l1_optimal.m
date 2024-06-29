function sol = l1_optimal(prob_par,optcost,wt_trgt)
    global nx nu Ad Bd d u_min u_max stage_cost alf_tp
    n = prob_par.n;
    N = prob_par.N;
    z_final = prob_par.zf;
    z_init = prob_par.z0;
    % Nmin = min(N);
    Nmax = max(N);
    cvx_begin
        variables x(nx,Nmax,n) u(nu,Nmax-1,n)
        expressions obj_val(n) l1_cost
        l1_cost = 0;
        for i=1:n
            obj_val(i) = 0;
            nn = 1:n; nn(i) = [];
            for j=1:Nmax-1
                if j<=N(i)-1
                    x(:,j+1,i) == Ad*x(:,j,i) + Bd*u(:,j,i) + d;
                    norm(u(:,j,i)) <= u_max;
                    u(3,j,i) >= u_min;
                    norm(u(1:2,j,i)) <= alf_tp*u(3,j,i);
                    obj_val(i) = obj_val(i) + stage_cost(x(:,j+1,i),u(:,j,i));
                    
                    for l=1:n-1
                        if j<=min(N(i),N(nn(l)))
                            l1_cost = l1_cost + norm(diag(wt_trgt{i,nn(l)}(:,j))*(x(:,j+1,i) - x(:,j+1,nn(l))),1);
                        end
                    end
                    
                else
                   x(:,j+1,i) == 0;
                   u(:,j,i) == 0;
                end
            end
            obj_val(i) <= (1+prob_par.opt_tol)*optcost(i)
            x(:,1,i) == z_init;
            x(:,N(i),i) == z_final(:,i);
        end
        minimize l1_cost
    cvx_end
    sol = struct;
    sol.X = cell(1,n);
    sol.U = cell(1,n);
    sol.cost = zeros(1,n);
    sol.l1cost = cvx_optval;
    for i=1:n
        sol.X{i} = x(:,1:N(i),i);
        sol.U{i} = u(:,1:N(i)-1,i);
        sol.cost(:) = obj_val(:);
    end
end