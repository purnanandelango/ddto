clearvars
close all
clc

zinit = [0;0;10;0;0;0;0;0];
rvtarg = [-15,  5, -20, -25;
           28, 28,  15,   0; 
           10,  5,  10,  10;
            1,  0,   0,   0;
            0,  1,   0,   1;
            0,  0,   0,   0];
cost_bound = 1000*[1,1,1,1];

prb = problem_data(04, ...          % K
                   030, ...         % scp_iters
                   5e1, ...         % w_ep
                   1.00, ...        % w_px
                   0.01, ...        % cost_factor
                   ...
                   zinit, ...
                   rvtarg, ...
                   cost_bound ...
                   );


load('recent_solution','xbar','ubar','taubar');
[xbar,ubar] = misc.create_initialization(prb,2, ...
                                         xbar,ubar,taubar);


[xbar,ubar] = scp.ctscvx_noparam(xbar,ubar,prb,@sys_cnstr_cost);
% [xbar,ubar] = scp.ctscvx_dvar_noparam(xbar,ubar,prb,@sys_cnstr_cost);

taubar = prb.tau;

% Simulate solution on fine grid
[tau,x,u] = disc.simulate_dyn(xbar(:,1),{prb.tau,ubar},@(t,x,u) prb.dyn_func(t,x,u),[0,1],prb.Kfine,prb.disc,prb.ode_solver);

[rbar,vbar,pbar,ybar,Tbar,sbar] = disassemble(xbar,ubar,prb.n,prb.ntarg+1);
[r,v,p,y,T,s] = disassemble(x,u,prb.n,prb.ntarg+1);

tvecbar(prb.K,prb.ntarg+1) = 0;
tvec(prb.Kfine,prb.ntarg+1) = 0;
for j = 1:prb.ntarg+1
    tvecbar(:,j) = disc.time_grid(prb.disc,prb.tau,sbar(:,j));
    tvec(:,j) = disc.time_grid(prb.disc,tau,s(:,j));
end

save('recent_solution','r','v','p','y','T','s','x','u','tvec','tau', ...
                       'rbar', 'vbar', 'pbar', 'ybar', 'Tbar', 'sbar', 'xbar', 'ubar', 'tvecbar', 'taubar', ...
                       'prb');

% plot_solution;