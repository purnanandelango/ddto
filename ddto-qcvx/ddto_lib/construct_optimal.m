function solopt = construct_optimal(prob_par)
    n = prob_par.n;
    N = prob_par.N;
    zf = prob_par.zf;
    z0 = prob_par.z0;
    
    fprintf("\nOptimal Solutions:\n")
    Xopt = cell(1,n);
    Uopt = cell(1,n);
    soln_cost = zeros(1,n);
    for i = 1:n
        [Xopt{i},Uopt{i},solve_status,soln_cost(i)] = compute_optimal(z0,zf(:,i),N(i));
        fprintf(" Target: %d, Cost: %.3f, Solve status: %s\n",prob_par.tag(i),soln_cost(i),solve_status);
        if ~strcmp("Solved",solve_status)
            error(horzcat('Target ',num2str(i),' is not reachable.'));
        end
    end
    solopt = struct;
    solopt.X = Xopt;
    solopt.U = Uopt;
    solopt.cost = soln_cost;
end

function [x,u,cvx_status,cvx_optval] = compute_optimal(z_init,z_final,N)
global nx nu Ad Bd d u_min u_max stage_cost alf_tp
    cvx_begin quiet
        variables x(nx,N) u(nu,N-1)
        expression obj_val(N-1)
        for j=1:N-1
            obj_val(j) = stage_cost(x(:,j+1),u(:,j));
            x(:,j+1) == Ad*x(:,j) + Bd*u(:,j) + d;
            norm(u(:,j)) <= u_max;
            u(3,j) >= u_min;
            norm(u(1:2,j)) <= alf_tp*u(3,j);
        end
        x(:,end) == z_final;
        x(:,1) == z_init;
    minimize sum(obj_val)
    cvx_end
end