function [] = plot_solution(prb,sol)

    interpreter = "tex";
    % interpreter = "latex";

    figure
    hold on
    idx_vec = []; 
    for j = setdiff(1:prb.n,prb.i)
        idx_vec(end+1) = sol.idx{j};
        plot3(sol.X{j}(1,sol.idx{j}:end),sol.X{j}(2,sol.idx{j}:end),sol.X{j}(3,sol.idx{j}:end),'.-','Color',prb.colors{j});%,'LineWidth',1.5,'MarkerSize',15)
    end
    plot3(sol.X{prb.i}(1,sol.idxmax:end),sol.X{prb.i}(2,sol.idxmax:end),sol.X{prb.i}(3,sol.idxmax:end),'.-','Color',prb.colors{prb.i});%,'LineWidth',1.5,'MarkerSize',15);

    idx_vec = sort(idx_vec);
    for j = 1:length(idx_vec)
        if j==1
            plot3(sol.X{prb.i}(1,1:idx_vec(j)),sol.X{prb.i}(2,1:idx_vec(j)),sol.X{prb.i}(3,1:idx_vec(j)),...
              '-','Color',prb.colors_trunk{j},'LineWidth',1.5,'MarkerSize',7,'Marker',prb.marker{j});                
        else
            plot3(sol.X{prb.i}(1,idx_vec(j-1):idx_vec(j)),sol.X{prb.i}(2,idx_vec(j-1):idx_vec(j)),sol.X{prb.i}(3,idx_vec(j-1):idx_vec(j)),...
              '-','Color',prb.colors_trunk{j},'LineWidth',1.5,'MarkerSize',7,'Marker',prb.marker{j});        
        end
    end

    if interpreter == "tex"
        annotation("textbox",'Position',[0.616,0.306142857142861,0.088571428571429,0.076428571428571],'EdgeColor','none','String','{\itz}^4');
        annotation("textbox",'Position',[0.647,0.531142857142861,0.088571428571429,0.076428571428571],'EdgeColor','none','String','{\itz}^3');
        annotation("textbox",'Position',[0.464000000000002,0.184476190476192,0.088571428571428,0.076428571428571],'EdgeColor','none','String','{\itz}^2');
        annotation("textbox",'Position',[0.362000000000001,0.163047619047621,0.088571428571429,0.076428571428571],'EdgeColor','none','String','{\itz}^1');
        annotation("textbox",'Position',[0.327428571428572,0.809952380952385,0.088571428571429,0.076428571428571],'EdgeColor','none','String','{\itz}^0');
    elseif interpreter == "latex"
        annotation("textbox",'Position',[0.678571428571428,0.296857142857141,0.088571428571429,0.076428571428571],'EdgeColor','none','String','$z^4$');
        annotation("textbox",'Position',[0.681428571428571,0.366857142857143,0.088571428571429,0.076428571428571],'EdgeColor','none','String','$z^3$');
        annotation("textbox",'Position',[0.468571428571429,0.255428571428572,0.088571428571429,0.076428571428571],'EdgeColor','none','String','$z^2$');
        annotation("textbox",'Position',[0.318571428571429,0.245428571428572,0.088571428571429,0.076428571428571],'EdgeColor','none','String','$z^1$');
        annotation("textbox",'Position',[0.201428571428571,0.714,0.088571428571429,0.076428571428571],'EdgeColor','none','String','$z^0$');
    end    

    grid on
    ax = gca;
    ax.PlotBoxAspectRatio = [1,1,1];
    ax.DataAspectRatio = [1,1,1];
    ax.View = [75,32];

    ax.ZLim = [-3,33];
    ax.XLim = [-1,41];
    ax.YLim = [-7,41];
    ax.XLabel.String = '[m]';
    ax.XLabel.Position = [20.78595919863642,-7.780066012867167,-7.998080814999412];
    ax.YLabel.String = '[m]';
    ax.YLabel.Position = [41.95130426648234,15.663663305361013,-8.327915792133012];
    ax.ZLabel.String = '[m]';      
end

