clearvars
close all
clc

rvtarg = [ 10, -10, -30, -15;
           30,  35,  15, -15; 
           10,  10,  10,  10;
            1,   0,   0,   0;
            0,   1,   0,   1;
            0,   0,   0,   0];
cost_bound = 1100*[1,1,1,1];

% M = 12;
% zinit = [ 10;
%          -10;
%           10;
%            0;
%            0;
%            0;
%            0;
%            0];
% 
% prb = problem_data_1(M, ...           % K
%                      020, ...         % scp_iters
%                      5e1, ...         % w_ep
%                      1.00, ...        % w_px
%                      0.001, ...        % cost_factor
%                      ...
%                      zinit, ...
%                      rvtarg, ...
%                      cost_bound ...
%                      );
% 
% load('results/trunk_1.mat','xbar','ubar','taubar');
% [xbar,ubar] = misc.create_initialization(prb,2, ...
%                                          xbar,ubar,taubar);

% load('results/trunk_1.mat','branch_point');
% M = 06;
% zinit = branch_point;
% rvtarg(:,end) = [];
% 
% prb = problem_data_2(M, ...           % K
%                      020, ...         % scp_iters
%                      5e1, ...         % w_ep
%                      1.00, ...        % w_px
%                      0.01, ...        % cost_factor
%                      ...
%                      zinit, ...
%                      rvtarg, ...
%                      cost_bound ...
%                      );
% 
% load('results/trunk_2.mat','xbar','ubar','taubar');
% [xbar,ubar] = misc.create_initialization(prb,2, ...
%                                          xbar,ubar,taubar);

load('results/trunk_2.mat','branch_point');
M = 03;
zinit = branch_point;
rvtarg(:,end-1:end) = [];

prb = problem_data_3(M, ...           % K
                     020, ...         % scp_iters
                     5e1, ...         % w_ep
                     1.00, ...        % w_px
                     0.02, ...        % cost_factor
                     ...
                     zinit, ...
                     rvtarg, ...
                     cost_bound ...
                     );

load('results/trunk_3.mat','xbar','ubar','taubar');
[xbar,ubar] = misc.create_initialization(prb,2, ...
                                         xbar,ubar,taubar);

% load('recent_solution','xbar','ubar','taubar');
% [xbar,ubar] = misc.create_initialization(prb,2, ...
%                                          xbar,ubar,taubar);

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

branch_point = xbar(1:2*prb.n+2,end);

save('recent_solution','r','v','p','y','T','s','x','u','tvec','tau', ...
                       'rbar', 'vbar', 'pbar', 'ybar', 'Tbar', 'sbar', 'xbar', 'ubar', 'tvecbar', 'taubar', ...
                       'prb','branch_point');

% plot_solution;