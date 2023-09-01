function [idx_x,idx_u,idx_s,idx_bet] = disassemble_traj_idx(n_states,n_ctrl,n_cnstr,n_targ)
    idx_x = zeros(n_states,n_targ);
    idx_u = zeros(n_ctrl,n_targ);
    idx_s = zeros(1,n_targ);
    idx_bet = [];
    if n_cnstr == 0
        for j = 1:n_targ
            idx_x(:,j) = (j-1)*n_states+1 : j*n_states;    
        end
    else
        idx_bet = zeros(n_cnstr,n_targ);
        for j = 1:n_targ
            idx_x(:,j)   = (j-1)*(n_states+n_cnstr)+1          : (j-1)*(n_states+n_cnstr)+n_states;
            idx_bet(:,j) = (j-1)*(n_states+n_cnstr)+n_states+1 : j*(n_states+n_cnstr);
        end 
    end
    for j = 1:n_targ
        idx_u(:,j) = (j-1)*(n_ctrl+1)+1 : (j-1)*(n_ctrl+1)+n_ctrl;
        idx_s(j)   = j*(n_ctrl+1);    
    end
end