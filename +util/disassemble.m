function [x,u,s,t,cum_traj_cost] = disassemble(tau,z,nu,n_cnstr,n_targ, ...
                                               time_grid_func,stage_cost_func)
% z  : vector containing states and constraint violation for all trajectories arranged as column vectors
% nu : vector containing control and dilation for all trajectories arranged as column vectors

    N               = size(z,2);
    n_states        = size(z,1)/n_targ - n_cnstr;
    n_ctrl          = size(nu,1)/n_targ - 1; 

    x               = zeros(n_states,N,n_targ);     % State
    u               = zeros(n_ctrl,N,n_targ);       % Control input
    s               = zeros(N,n_targ);              % Dilation factor
    t               = zeros(N,n_targ);              % Time grid
    cum_traj_cost   = zeros(N,n_targ);              % Cumulative trajectory cost

    [idx_x,idx_u,idx_s] = util.disassemble_traj_idx(n_states,n_ctrl,n_cnstr,n_targ);
    for j = 1:n_targ
        for k = 1:N
            x(:,k,j) = z(idx_x(:,j),k);
            u(:,k,j) = nu(idx_u(:,j),k);
            s(k,j) = nu(idx_s(j),k);
            if k < N
                cum_traj_cost(k+1,j) = cum_traj_cost(k,j) + stage_cost_func(x(:,k,j),u(:,k,j));    
            end
        end
        t(:,j) = time_grid_func(tau,x(:,:,j),s(:,j));
    end

end