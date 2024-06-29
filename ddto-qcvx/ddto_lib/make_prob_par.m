function prob_par = make_prob_par(z0,zf,n,N,tag,opt_tol)
    prob_par = struct;
    prob_par.z0 = z0;                % initial state
    prob_par.zf = zf;                % target states   
    prob_par.N = N;                  % no. of discretization points 
    prob_par.n = n;                  % no. of targets   
    prob_par.tag = tag;              % target tag
    prob_par.opt_tol = opt_tol;      % optimality tolerance
end