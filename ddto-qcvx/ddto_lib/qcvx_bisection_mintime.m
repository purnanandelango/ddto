function sol = qcvx_bisection_mintime(prob_par)

    n = prob_par.n;
    N = prob_par.N;
    
    sol = struct;
    sol.X = cell(1,n);
    sol.U = cell(1,n);
    sol.cost = zeros(1,n);
    sol.N = zeros(1,n);

    fprintf("Computation of minimum-time solutions\n-------------------------------------\n")
    for i = 1:n
        tau_max = N(i);
        tau_min = 1;

        fprintf("\nBisection Search:\n Iteration: %.2d,  Nmin: %.2d, Nmax: %.2d, Status: %s",0,tau_min,tau_max,"None");

        % bisection search for solving quasiconvex optimization problem
        iter_count = 1;
        while tau_max-tau_min>1 % termination criteria   

           tau = ceil(0.5*(tau_max+tau_min)); % update tau

           [~,~,solve_status,~] = compute_feasible_mintime(prob_par,i,tau);

           if strcmp(string(solve_status),"Solved") % problem is feasible
                tau_max = tau;
           else
                tau_min = tau;    
           end 

           fprintf("\n Iteration: %.2d,  Nmin: %.2d, Nmax: %.2d, Status: %s",iter_count,tau_min,tau_max,solve_status);

           iter_count = iter_count + 1;

        end
        sol.N(i) = tau_max;
        fprintf("\n Terminated (Nmax-Nmin=1).")
        [sol.X{i},sol.U{i},solve_status,sol.cost(i)] = compute_feasible_mintime(prob_par,i,sol.N(i));
        if strcmp(solve_status,"Solved")
            fprintf(" Nstr = %d",sol.N(i))        
        else
            fprintf("\n")
            error("Bisection search unsuccessful. Problem is unsolved.")
        end
    end
    
    fprintf("\n\nMinimum-time Solutions:\n")
    for i=1:n
        fprintf(" Target: %d, N: %2d, Nstr: %2d, Cost: %.3f\n",prob_par.tag(i),prob_par.N(i),sol.N(i),sol.cost(i));
    end

end
function [x,u,cvx_status,soln_cost] = compute_feasible_mintime(prob_par,i,k)
    global Ad Bd d nx nu stage_cost u_min u_max alf_tp
    z_init = prob_par.z0;
    z_final = prob_par.zf(:,i);
    N = prob_par.N(i);

    cvx_begin quiet
        variables x(nx,N) u(nu,N-1)
        find [x,u];
        expression obj_val(N-1)
        for j=1:N-1
            x(:,j+1) == Ad*x(:,j) + Bd*u(:,j) + d;
            norm(u(:,j)) <= u_max;
            u(3,j) >= u_min;
            norm(u(1:2,j)) <= alf_tp*u(3,j);
            obj_val(j) = stage_cost(x(:,j+1),u(:,j));
        end
        for l=k:N
            x(:,l) == z_final;
        end
        x(:,1) == z_init;
    cvx_end
    soln_cost = sum(obj_val(1:k));
    
end