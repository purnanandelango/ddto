clearvars
close all
clc

prb = problem_data(1,1000);

sol  = solve_ddto_micp(prb);
sol1 = solve_ddto_micp_relax_v1(prb,[]);
sol2 = solve_ddto_micp_relax_v2(prb,1,[]);

plot_solution(prb,sol);
plot_solution(prb,sol1);
plot_solution(prb,sol2);
