clearvars
clc

prb = problem_data_2D(14, ...           % K
                      10, ...           % scp_iters
                      5e1, ...          % wvc
                      5e1, ...          % wvb    
                      0.01, ...         % wtr
                      0.05);            % cost_factor

load('recent_solution','xbar','ubar');
[xbar,ubar] = misc.create_initialization(prb,1, ...
                                         xbar,ubar);

[xbar,ubar] = scp.run_ptr_noparam(xbar,ubar,prb,@sys_cnstr_cost);

% Simulate solution on fine grid in [0,1]
[tau,x,u] = disc.simulate_dyn(xbar(:,1),{prb.tau,ubar},@(t,x,u) prb.dyn_func(t,x,u),[0,1],prb.Kfine,prb.disc);

[r,   v,   T,   s,   tvec,   nrm_v,   nrm_T,...
 rbar,vbar,Tbar,sbar,tvecbar,nrm_vbar,nrm_Tbar,...
 traj_cost]                                         = util.disassemble(x,u,tau,xbar,ubar,prb);    

save('recent_solution','r','v','nrm_v','T','nrm_T','s','tvec','x','u','prb','tau',...
                       'rbar','vbar','nrm_vbar','Tbar','nrm_Tbar','sbar','tvecbar','xbar','ubar','traj_cost');

plot_solution;