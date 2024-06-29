function soldd = qcvx_bisection_dd(prob_par,optcost)

    n = prob_par.n;
    N = prob_par.N;

    tau_max = min(N)-2;
    tau_min = 0;
    
    fprintf("\nBisection Search:\n Iteration: %.2d,  Nmin: %.2d, Nmax: %.2d, Status: %s",0,tau_min,tau_max,"None");
    
    % bisection search for solving quasiconvex optimization problem
    iter_count = 1;
    while tau_max-tau_min>1 % termination criteria   
        
       tau = ceil(0.5*(tau_max+tau_min)); % update tau
       
       [~,~,solve_status,~,~] = compute_feasible_ddto(prob_par,optcost,tau,false);
       
       if strcmp(string(solve_status),"Solved") % problem is feasible
            tau_min = tau;
       else
            tau_max = tau;    
       end 
       
       fprintf("\n Iteration: %.2d,  Nmin: %.2d, Nmax: %.2d, Status: %s",iter_count,tau_min,tau_max,solve_status);

       iter_count = iter_count + 1;
       
    end
    tau_str = tau_min;
    fprintf("\n Terminated (Nmax-Nmin=1).")
    [X,U,solve_status,soln_cost,soln_cost_dd] = compute_feasible_ddto(prob_par,optcost,tau_str,false);
    if strcmp(solve_status,"Solved")
        fprintf(" Nstr = %d",tau_str)        
    else
        fprintf("\n")
        error("Bisection search unsuccessful. Problem is unsolved.")
    end
    fprintf("\n\nDeferred-decision Solutions:\n")
    for i=1:n
        fprintf(" Target: %d, Cost: %.3f\n",prob_par.tag(i),soln_cost(i));
    end
    soldd = struct;
    soldd.X = X;
    soldd.U = U;
    soldd.cost = soln_cost;
    soldd.cost_dd = soln_cost_dd; 
    soldd.dd_idx = tau_str;
end

function [X,U,slv_status,soln_cost,soln_cost_dd] = compute_feasible_ddto(prob_par,optcost,k,print_flag)
global nx nu Ad Bd d u_min u_max stage_cost alf_tp

    n = prob_par.n;
    N = prob_par.N;
    eps_tol = prob_par.opt_tol;
    z_init = prob_par.z0;
    z_final = prob_par.zf;
    
    Nmax = max(N);

%     yalmip clear
%     x = sdpvar(nx,Nmax,n);
%     u = sdpvar(nu,Nmax-1,n);
%     subopt = sdpvar(Nmax-1,n);
%     cnstr = [];
%     for i = 1:n
%         nn = 1:n; nn(i) = [];
%         for j = 1:Nmax-1
%             if j<=N(i)-1
%                 cnstr = [cnstr;
%                          x(:,j+1,i) == Ad*x(:,j,i) + Bd*u(:,j,i) + d;
%                          norm(u(:,j,i)) <= u_max;
%                          u(3,j,i) >= u_min;
%                          norm(u(1:2,j,i)) <= alf_tp*u(3,j,i)];
% 
%                 %%% SUB-OPTIMALITY
%                 subopt(j,i) = stage_cost(x(:,j+1,i),u(:,j,i));  
% 
%                 %%% IDENTICALITY 
%                 if k>0
%                     if j<=k 
%                         for l = 1:n-1
%                             cnstr = [cnstr;
%                                      u(:,j,i) == u(:,j,nn(l))]; 
%                         end
%                     end
%                 end
% 
%             else
%                 cnstr = [cnstr;
%                          subopt(j,i) == 0;
%                          x(:,j+1,i) == zeros(nx,1);
%                          u(:,j,i) == zeros(nu,1)];                
%             end
%         end
%         %%% SUB-OPTIMALITY          
%         cnstr = [cnstr;
%                  sum(subopt(:,i)) + optcost.dd <= (1+eps_tol)*optcost.opt(i);
%                  x(:,N(i),i) == z_final(:,i);
%                  x(:,1,i) == z_init];
%     end
%     yal_out = optimize(cnstr,0,sdpsettings('solver','mosek'));
%     x = value(x);
%     u = value(u);
%     % subopt = value(subopt);
%     if print_flag
%         fprintf("\n Solve status : %s\n",yalmiperror(yal_out.problem))
%     end
%     if yal_out.problem == 0
%         slv_status = "Solved";
%     else
%         slv_status = "Unsolved";
%     end
    
    cvx_begin quiet
        variables x(nx,Nmax,n) u(nu,Nmax-1,n)
        find [x,u];
        expression subopt(Nmax-1,n);
        for i=1:n
            nn = 1:n; nn(i) = [];
            for j=1:Nmax-1
                if j<=N(i)-1
                   x(:,j+1,i) == Ad*x(:,j,i) + Bd*u(:,j,i) + d;
                   norm(u(:,j,i)) <= u_max;
                   u(3,j,i) >= u_min;
                   norm(u(1:2,j,i)) <= alf_tp*u(3,j,i);
                   
                   %%% SUB-OPTIMALITY
                   subopt(j,i) = stage_cost(x(:,j+1,i),u(:,j,i));                                                      

                   %%% IDENTICALITY 
                   if k>0
                       if j<=k 
                           for l = 1:n-1
                                u(:,j,i) == u(:,j,nn(l)); 
                           end
                       end
                   end
                   
                else
                    subopt(j,i) = 0;
                    x(:,j+1,i) == zeros(nx,1);
                    u(:,j,i) == zeros(nu,1);
                end
            end
            
            %%% SUB-OPTIMALITY          
            sum(subopt(:,i)) + optcost.dd <= (1+eps_tol)*optcost.opt(i);
            
            x(:,N(i),i) == z_final(:,i);
            x(:,1,i) == z_init;            
        end
    cvx_end
    if print_flag
        fprintf("\n Solve status : %s\n",cvx_status)
    end
    slv_status  = cvx_status;  

    X = cell(1,n); 
    U = cell(1,n); 
    soln_cost = zeros(1,n);
    soln_cost_dd = [];
    for i = 1:n
        X{i} = x(:,1:N(i),i);
        U{i} = u(:,1:N(i)-1,i);
        for j=1:N(i)-1
            soln_cost(i) = soln_cost(i) + stage_cost(x(:,j+1,i)-z_final(:,i),u(:,j,i));
            if j==k && i==1
                soln_cost_dd(end+1) = soln_cost(i);
            end
        end
    end
end