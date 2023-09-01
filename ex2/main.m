clearvars
clc

prb = problem_data_2D(14, ...           % K
                      10, ...           % scp_iters
                      5e1, ...          % wvc
                      0.01, ...         % wtr
                      0.10);            % cost_factor

load('recent_solution','zbar','nubar');
[zbar,nubar] = misc.create_initialization(prb,2, ...
                                          zbar,nubar,[]);

[zbar,nubar] = scp.run_ptr_noparam(zbar,nubar,prb,@sys_cnstr_cost);
taubar = prb.tau;

% Simulate solution on fine grid in [0,1]
[tau,z,nu] = disc.simulate_dyn(zbar(:,1),{prb.tau,nubar},@(t,z,nu) prb.dyn_func(t,z,nu),[0,1],prb.Kfine,prb.disc);

[x,u,s,tvec] = util.disassemble(tau,z,nu,prb.m,prb.ntarg,prb.time_grid,prb.cost_term);
[xbar,ubar,sbar,tvecbar,cum_traj_cost] = util.disassemble(taubar,zbar,nubar,prb.m,prb.ntarg,prb.time_grid,prb.cost_term);

nrmv(prb.Kfine,prb.ntarg) = 0;
nrmu(prb.Kfine,prb.ntarg) = 0;
nrmvbar(prb.K,prb.ntarg) = 0;
nrmubar(prb.K,prb.ntarg) = 0;
r    = x(1:prb.n,:,:);
v    = x(prb.n+1:2*prb.n,:,:);
rbar = xbar(1:prb.n,:,:);
vbar = xbar(prb.n+1:2*prb.n,:,:);
for j = 1:prb.ntarg
    nrmv(:,j) = misc.column_func(v(:,:,j),@norm);
    nrmu(:,j) = misc.column_func(u(:,:,j),@norm);
    nrmvbar(:,j) = misc.column_func(vbar(:,:,j),@norm);
    nrmubar(:,j) = misc.column_func(ubar(:,:,j),@norm);  
    fprintf('\nTarget %d - Final position error: %.3f\n         - Final velocity error: %.3f\n         -      Trajectory cost: %.3f\n',j,norm(r(:,end,j)-prb.rK(:,j)),norm(v(:,end,j)-prb.vK(:,j)),cum_traj_cost(end,j));    
end

save('recent_solution','r','v','nrmv','u','nrmu','s','tvec','z','nu','tau','prb',...
                       'rbar','vbar','nrmvbar','ubar','nrmubar','sbar','tvecbar','zbar','nubar','taubar','cum_traj_cost');

plot_solution;