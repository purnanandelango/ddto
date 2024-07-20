clearvars
close all
clc

prb = problem_data(1,100);

sol = solve_ddto_mip(prb);

colors = {[1,0,0.5],[0,1,0.5],[0.5,0,1],[0.5,0.5,1]};

figure
hold on
for j = setdiff(1:prb.n,prb.i)
    plot3(sol.X{j}(1,:),sol.X{j}(2,:),sol.X{j}(3,:),'.-','Color',colors{j});
end
plot3(sol.X{prb.i}(1,:),sol.X{prb.i}(2,:),sol.X{prb.i}(3,:),'.-','Color',[0,0,0]);
grid on
ax = gca;
ax.PlotBoxAspectRatio = [1,1,1];
ax.DataAspectRatio = [1,1,1];
ax.View = [186,20];