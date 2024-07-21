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

weights = [];
for k = 1:1
    sol2 = solve_ddto_micp_relax_v2(prb,weights,[]);
    weights = sol2.Xi;
end
sol2.objval
 
plot_solution(prb,sol);
plot_solution(prb,sol2);
plot_solution(prb,sol1);