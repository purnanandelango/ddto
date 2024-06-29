function plot_soln(soldd,solopt,prob_par,plt_mrk,mrk_target)
    
    n = prob_par.n;
    ax_lim = [-0.5*max(max(prob_par.zf([1,2,3],:))),max(max(prob_par.zf([1,2,3],:)))];
    
    for j=1:n
    % primary trajectory
    plot3(soldd.X{j}(1,:),soldd.X{j}(2,:),soldd.X{j}(3,:),'LineStyle',plt_mrk{1}{1},'Color',plt_mrk{1}{2},'Marker',plt_mrk{1}{3},'LineWidth',plt_mrk{1}{4},'DisplayName','DD Traj.');
    if j == 1
        hold on
    end
    if isstruct(solopt)
        plot3(solopt.X{j}(1,:),solopt.X{j}(2,:),solopt.X{j}(3,:),'LineStyle',plt_mrk{2}{1},'Color',plt_mrk{2}{2},'Marker',plt_mrk{2}{3},'LineWidth',plt_mrk{2}{4},'DisplayName','Optimal Traj.');
    end
    legend('AutoUpdate','off');
    end
    legend('AutoUpdate','on','location','west');

    if mrk_target
        m_clr = {[0,0.5,0.5],[0.5,0.5,0],[0.5,0,0.5],[0,0.5,0],[0.5,0.5,0.5],[0.2,0.4,0]};
        for i=1:n
        plot(prob_par.zf(1,i),prob_par.zf(2,i),'Marker','o','Color',m_clr{i},'DisplayName',horzcat('Target ',num2str(prob_par.tag(i))));    
        end
        plot(prob_par.z0(1),prob_par.z0(2),'ok','DisplayName','Source');
    end
    xlim(ax_lim)
    ylim(ax_lim)
    axis square
end