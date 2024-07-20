function [] = plot_solution(prb,sol)
    figure
    hold on
    plot3(sol.X{prb.i}(1,:),sol.X{prb.i}(2,:),sol.X{prb.i}(3,:),'.-','Color',[0,0,0],'LineWidth',3,'MarkerSize',20);
    for j = setdiff(1:prb.n,prb.i)
        plot3(sol.X{j}(1,:),sol.X{j}(2,:),sol.X{j}(3,:),'.-','Color',prb.colors{j},'LineWidth',1,'MarkerSize',15);
    end
    grid on
    ax = gca;
    ax.PlotBoxAspectRatio = [1,1,1];
    ax.DataAspectRatio = [1,1,1];
    ax.View = [76,34];
end

