clear all
close all
clc 

addpath('./ddto_lib')
addpath('./plot_lib')

% Problem parameters
set_params;

% Specification of targets
target_order = [2,3,1,4]; % Order of target rejection
target_set1 = prob_par.tag;
target_set2 = prob_par.tag; target_set2(target_order(1)) = [];
target_set3 = prob_par.tag; target_set3(target_order([1,2])) = [];

%% Determination of minimum-time to reach all targets
% solmint = qcvx_bisection_mintime(prob_par);
% prob_par.N = min(ceil(solmint.N*1.5),max(prob_par.N));

prob_par.N = 20*[1,1,1,1];
% prob_par.N = [21,18,19,20];
% prob_par.N = [12,11,10,12];

%% Optimal solution from source
optcost = struct;
solopt = construct_optimal(prob_par); 
optcost.opt = solopt.cost;

%% First branching
optcost.dd = 0;
soldd = qcvx_bisection_dd(prob_par,optcost);

%% Second branching
% Setup
prob_par2 = make_prob_par(soldd.X{1}(:,soldd.dd_idx+1),...              % z0
                          prob_par.zf(:,target_set2),...                % zf
                          prob_par.n-1,...                              % n
                          prob_par.N(target_set2)-soldd.dd_idx,...      % N    
                          target_set2,...                               % Target tag    
                          1.0);
                          % prob_par.opt_tol);                            % Optimality tolerance                      

optcost.dd = soldd.cost_dd;
soldd2 = qcvx_bisection_dd(prob_par2,optcost);

%% Third branching
% Setup
prob_par3 = make_prob_par(soldd2.X{1}(:,soldd2.dd_idx+1),...                        % z0
                          prob_par.zf(:,target_set3),...                            % zf
                          prob_par2.n-1,...                                         % n
                          prob_par.N(target_set3)-soldd2.dd_idx-soldd.dd_idx,...    % N
                          target_set3,...                                           % Target tag
                          1.0);
                          % prob_par.opt_tol);                                        % Optimality tolerance

optcost.dd = soldd.cost_dd+soldd2.cost_dd;
soldd3 = qcvx_bisection_dd(prob_par3,optcost);

%% l1 optimal solutions

% weight matrix
% [11 12 13 14
%  21 22 23 24
%  31 32 33 34
%  41 42 43 44]

% wt_mat = [0,    1,  100,    1000;
%           1,    0,  1,      1;
%           100,  1,  0,      100;
%           1000, 1,  100,    0];
% 
% wt_trgt = create_l1_wt(wt_mat,prob_par.N);
% 
% soll1 = l1_optimal(prob_par,optcost.opt,wt_trgt);

%% Re-weighted l1 optimal solution

% wt_trgt = create_l1_wt(ones(prob_par.n)-eye(prob_par.n),prob_par.N);
% wt_trgt = create_l1_wt(wt_mat,prob_par.N);
% for iter = 1:5
%     solre_l1 = l1_optimal(prob_par,optcost.opt,wt_trgt);
%     wt_trgt = update_l1_wt(solre_l1.X,wt_trgt,prob_par.N,1e-3); 
% end    

%% Plotting

plot_condensed(prob_par,soldd,prob_par2,soldd2,prob_par3,soldd3,target_order);
% plot_condensed_l1(soll1,{'.-','l1',[0.4,0,0.9]})
% plot_condensed_l1(solre_l1,{'.-','re-weight l1',[0.9,0.4,0.0]})

%% 
% figure
% plot_soln(soldd,solopt,prob_par,{{'-',[0,0,1],'.',1},{'-',[1,0,0],'.',1}},true);
% plot_soln(soldd2,NaN,prob_par2,{{'-',[0.5,0,1],'.',1.5}},false);
% plot_soln(soldd3,NaN,prob_par3,{{'-',[0,0.5,1],'.',2}},false);
% title(horzcat( num2str(100*prob_par.opt_tol) ,' \%, ', num2str(100*prob_par2.opt_tol) ,' \%, ', num2str(100*prob_par3.opt_tol) ,' \%'));

save("solution","optcost","solopt","prob_par","soldd","prob_par2","soldd2","prob_par3","soldd3","target_order");

rmpath('./ddto_lib')
rmpath('./plot_lib')

%% 

function wt_new = update_l1_wt(X,wt,N,tol_val)
    n = size(wt,1);
    wt_new = wt;
    for i = 1:n
       nn = 1:n;
       nn(i) = [];
       for j = 1:max(N)-1
            for l = 1:n-1
                if j<=min([N(i),N(nn(l))])-1
                    wt_new{i,nn(l)}(:,j) = 1 ./ ( abs(X{i}(:,j+1)-X{nn(l)}(:,j+1)) + tol_val );
                end
            end
       end
    end
end

function wt_trgt = create_l1_wt(wt_mat,N)
% create the weights cell array for the l1 ddto
global nx
    n = size(wt_mat,1);
    wt_trgt = cell(n,n);
    for i=1:n
        for j=1:n
            if i==j
                wt_trgt{i,j} = [];
            else
                wt_trgt{i,j} = repmat(wt_mat(i,j)*ones(nx,1),[1,min(N([i,j]))]);
            end
        end
    end
end
