clearvars
close all
clc

prb = problem_data(4,1e4);

sol  = solve_ddto_micp(prb);
sol.objval

prb.M = 1e5;

weights = [];
for k = 1:1
    sol1 = solve_ddto_micp_relax_v1(prb,weights);
    weights = sol1.Xi;
end
sol1.objval

sol1.idx = sol.idx;
sol1.idxmax = sol.idxmax;

% weights = [];
% for k = 1:1
%     sol2 = solve_ddto_micp_relax_v2(prb,weights,[]);
%     weights = sol2.Xi;
%     % [weights{1},weights{3},weights{4}]
% end
% sol2.objval

% sol2.idx = sol.idx;
% sol2.idxmax = sol.idxmax;

plot_solution(prb,sol);
exportgraphics(gca,"ddto-micp.pdf",'ContentType','vector');
savefig(gcf,'ddto-micp.fig'); 
% plot_solution(prb,sol1);
% plot_solution(prb,sol2);
% exportgraphics(gca,"ddto-relax.pdf",'ContentType','vector');