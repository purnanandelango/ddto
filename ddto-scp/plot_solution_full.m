clearvars
close all


sol1 = load("results/trunk_1.mat");
sol2 = load("results/trunk_2.mat");
sol3 = load("results/trunk_3.mat");

prb = sol1.prb;

ntarg_p_1 = prb.ntarg+1;

interpreter = "tex";        plt.setfig("tex");
% interpreter = "latex";      plt.setfig("latex")
saveplot = true;

Colors_branch = {[1,0.5,0.5],[0,0.7,0.3],[0.5,0.5,1],[1,0.65,0.3]};
Colors_trunk = {[0.1,0.1,0.2],[0.3,0.3,0.4],[0.5,0.5,0.6]};
Color_bound = [161,240,242]/255;

f1 = figure('Position',[100,100,700,700]);

% Plot trunk and branch position trace
plot3(sol1.r(1,:,1),     sol1.r(2,:,1),     sol1.r(3,:,1),     '-','Color',Colors_trunk{1});
hold on 
plot3(sol1.rbar(1,:,1),  sol1.rbar(2,:,1),  sol1.rbar(3,:,1),  'o','Color',Colors_trunk{1},'MarkerSize',7);
plot3(sol1.r(1,:,end),   sol1.r(2,:,end),   sol1.r(3,:,end),   '-','Color',Colors_branch{ntarg_p_1-1});
plot3(sol1.rbar(1,:,end),sol1.rbar(2,:,end),sol1.rbar(3,:,end),'.','Color',Colors_branch{ntarg_p_1-1});

plot3(sol2.r(1,:,1),     sol2.r(2,:,1),     sol2.r(3,:,1),     '-','Color',Colors_trunk{2});
plot3(sol2.rbar(1,:,1),  sol2.rbar(2,:,1),  sol2.rbar(3,:,1),  's','Color',Colors_trunk{2},'MarkerSize',7);
plot3(sol2.r(1,:,end),   sol2.r(2,:,end),   sol2.r(3,:,end),   '-','Color',Colors_branch{ntarg_p_1-2});
plot3(sol2.rbar(1,:,end),sol2.rbar(2,:,end),sol2.rbar(3,:,end),'.','Color',Colors_branch{ntarg_p_1-2});

plot3(sol3.r(1,:,1),       sol3.r(2,:,1),       sol3.r(3,:,1),       '-','Color',Colors_trunk{3});
plot3(sol3.rbar(1,:,1),    sol3.rbar(2,:,1),    sol3.rbar(3,:,1),    '^','Color',Colors_trunk{3},'MarkerSize',7);
plot3(sol3.r(1,:,end),     sol3.r(2,:,end),     sol3.r(3,:,end),     '-','Color',Colors_branch{ntarg_p_1-3});
plot3(sol3.rbar(1,:,end),  sol3.rbar(2,:,end),  sol3.rbar(3,:,end),  '.','Color',Colors_branch{ntarg_p_1-3});
plot3(sol3.r(1,:,end-1),   sol3.r(2,:,end-1),   sol3.r(3,:,end-1),   '-','Color',Colors_branch{ntarg_p_1-4});
plot3(sol3.rbar(1,:,end-1),sol3.rbar(2,:,end-1),sol3.rbar(3,:,end-1),'.','Color',Colors_branch{ntarg_p_1-4});

grid on

ax = gca;
ax.Box = "off";
ax.View = [74,16];
ax.ZLim = [3,17];
ax.XLim = [-32,12];
ax.YLim = [-16,35];
ax.XLabel.String = '[m]';
ax.XLabel.Position = [-11.140111463662151,-23.218241859117768,-0.116603058341084];
ax.YLabel.String = '[m]';
ax.YLabel.Position = [12.349562166848443,5.590137687474993,-1.93024396695246];
ax.ZLabel.String = '[m]';

ax.DataAspectRatio = [1,1,1];
ax.PlotBoxAspectRatio = [1,1,1]; 

if interpreter == "latex"
    annotation("textbox",'Position',[0.181385714285714,0.5004,0.089,0.0657],'EdgeColor','none','String','$z^4$');
    annotation("textbox",'Position',[0.494285714285714,0.578571428571428,0.088958544049944,0.065714285714286],'EdgeColor','none','String','$z^3$');
    annotation("textbox",'Position',[0.811428571428571,0.533285714285714,0.088958544049944,0.065714285714286],'EdgeColor','none','String','$z^2$');
    annotation("textbox",'Position',[0.824285714285714,0.459,0.088958544049944,0.065714285714286],'EdgeColor','none','String','$z^1$');
elseif interpreter == "tex"
    annotation("textbox",'Position',[0.185671428571428,0.513257142857143,0.089,0.0657],'EdgeColor','none','String','{\itz}^4');
    annotation("textbox",'Position',[0.498571428571428,0.591428571428571,0.088958544049944,0.065714285714286],'EdgeColor','none','String','{\itz}^3');
    annotation("textbox",'Position',[0.812857142857143,0.543285714285714,0.088958544049944,0.065714285714286],'EdgeColor','none','String','{\itz}^2');
    annotation("textbox",'Position',[0.822857142857143,0.471857142857143,0.088958544049944,0.065714285714286],'EdgeColor','none','String','{\itz}^1');
end

[X1,Y1,Z1] = ellipsoid(0,0,0,1/prb.robs{1}(1),1/prb.robs{1}(2),1/prb.robs{1}(3),30);
X1 = X1 + sol1.prb.qobs(1,1);
Y1 = Y1 + sol1.prb.qobs(2,1);
Z1 = Z1 + sol1.prb.qobs(3,1);
ellip1 = surf(X1,Y1,Z1,...
    ...% 'EdgeColor',[1,0.5,0.2],'EdgeAlpha',0.5,'FaceColor',[1,1,0],'FaceAlpha',0.2);
       'EdgeColor',[0.3,0,0.3],'EdgeAlpha',0.2,'FaceColor',[0,0,0],'FaceAlpha',0.05);
rotate(ellip1,[0,0,1],0);

[X2,Y2,Z2] = ellipsoid(0,0,0,1/prb.robs{2}(1),1/prb.robs{2}(2),1/prb.robs{2}(3),30);
X2 = X2 + sol1.prb.qobs(1,2);
Y2 = Y2 + sol1.prb.qobs(2,2);
Z2 = Z2 + sol1.prb.qobs(3,2);
ellip2 = surf(X2,Y2,Z2,...
    ...% 'EdgeColor',[1,0.5,0.2],'EdgeAlpha',0.5,'FaceColor',[1,1,0],'FaceAlpha',0.2);
       'EdgeColor',[0.3,0,0.3],'EdgeAlpha',0.2,'FaceColor',[0,0,0],'FaceAlpha',0.05);
rotate(ellip2,[0,0,1],0);

if interpreter == "latex"
    ticklab = ax.XTickLabel;
    for j = 1:length(ticklab)
        ticklab{j} = horzcat('$',ticklab{j},'$');
    end
    ax.XTickLabel = ticklab;
    ticklab = ax.YTickLabel;
    for j = 1:length(ticklab)
        ticklab{j} = horzcat('$',ticklab{j},'$');
    end
    ax.YTickLabel = ticklab;
    ticklab = ax.ZTickLabel;
    for j = 1:length(ticklab)
        ticklab{j} = horzcat('$',ticklab{j},'$');
    end
    ax.ZTickLabel = ticklab;
elseif interpreter == "tex"
    ax.XTickLabel = strrep(ax.XTickLabel,'-',char(8722));
    ax.YTickLabel = strrep(ax.YTickLabel,'-',char(8722));
    ax.ZTickLabel = strrep(ax.ZTickLabel,'-',char(8722));
end

if saveplot
    exportgraphics(ax,"results/position_"+interpreter+".pdf",'ContentType','vector');
    savefig(f1,'results/position.fig');
end

ToFmax = max(sol1.ToFtrunk+sol2.ToFtrunk+sol3.ToFtrunk+sol3.tvecbar(end,2:end));

f2 = figure('Position',[1500,101,1400,700],'Units','normalized');

% Thrust pointing angle
subplot(3,2,1)
plot(linspace(0,ToFmax,100),prb.deltamax*ones(1,100),'-','Color',Color_bound,'LineWidth',4);
hold on 
plot(              sol1.tvec   (:,  1),sol1.point_angle   (:,  1),'-','Color',Colors_trunk{1});
plot(              sol1.tvecbar(:,  1),sol1.point_anglebar(:,  1),'o','Color',Colors_trunk{1},'MarkerSize',7);
plot(sol1.ToFtrunk+sol1.tvec   (:,end),sol1.point_angle   (:,end),'-','Color',Colors_branch{ntarg_p_1-1});
plot(sol1.ToFtrunk+sol1.tvecbar(:,end),sol1.point_anglebar(:,end),'.','Color',Colors_branch{ntarg_p_1-1});

plot(sol1.ToFtrunk+              sol2.tvec   (:,  1),sol2.point_angle   (:,  1),'-','Color',Colors_trunk{2});
plot(sol1.ToFtrunk+              sol2.tvecbar(:,  1),sol2.point_anglebar(:,  1),'s','Color',Colors_trunk{2},'MarkerSize',7);
plot(sol1.ToFtrunk+sol2.ToFtrunk+sol2.tvec   (:,end),sol2.point_angle   (:,end),'-','Color',Colors_branch{ntarg_p_1-2});
plot(sol1.ToFtrunk+sol2.ToFtrunk+sol2.tvecbar(:,end),sol2.point_anglebar(:,end),'.','Color',Colors_branch{ntarg_p_1-2});

plot(sol1.ToFtrunk+sol2.ToFtrunk+              sol3.tvec   (:,  1),  sol3.point_angle   (:,  1),  '-','Color',Colors_trunk{3});
plot(sol1.ToFtrunk+sol2.ToFtrunk+              sol3.tvecbar(:,  1),  sol3.point_anglebar(:,  1),  '^','Color',Colors_trunk{3},'MarkerSize',7);
plot(sol1.ToFtrunk+sol2.ToFtrunk+sol3.ToFtrunk+sol3.tvec   (:,end),  sol3.point_angle   (:,end),  '-','Color',Colors_branch{ntarg_p_1-3});
plot(sol1.ToFtrunk+sol2.ToFtrunk+sol3.ToFtrunk+sol3.tvecbar(:,end),  sol3.point_anglebar(:,end),  '.','Color',Colors_branch{ntarg_p_1-3});
plot(sol1.ToFtrunk+sol2.ToFtrunk+sol3.ToFtrunk+sol3.tvec   (:,end-1),sol3.point_angle   (:,end-1),'-','Color',Colors_branch{ntarg_p_1-4});
plot(sol1.ToFtrunk+sol2.ToFtrunk+sol3.ToFtrunk+sol3.tvecbar(:,end-1),sol3.point_anglebar(:,end-1),'.','Color',Colors_branch{ntarg_p_1-4});

% title('Thrust pointing angle [deg]')
ylabel('[deg]');
xlabel('[s]');
xlim([0,ToFmax])
ylim([0,1.1*prb.deltamax]);
ax = gca;
ax.Box = "off";
ax.XLabel.Position = [4.933679040598237,-10.536389208570618,-0.999999999999986];
if saveplot
    exportgraphics(ax,"results/pointangle_"+interpreter+".pdf");
end

subplot(3,2,2)
plot(linspace(0,ToFmax,100),prb.vmax*ones(1,100),'-','Color',Color_bound,'LineWidth',4);
hold on 
plot(              sol1.tvec   (:,  1),sol1.nrm_v   (:,  1),'-','Color',Colors_trunk{1});
plot(              sol1.tvecbar(:,  1),sol1.nrm_vbar(:,  1),'o','Color',Colors_trunk{1},'MarkerSize',7);
plot(sol1.ToFtrunk+sol1.tvec   (:,end),sol1.nrm_v   (:,end),'-','Color',Colors_branch{ntarg_p_1-1});
plot(sol1.ToFtrunk+sol1.tvecbar(:,end),sol1.nrm_vbar(:,end),'.','Color',Colors_branch{ntarg_p_1-1});

plot(sol1.ToFtrunk+              sol2.tvec   (:,  1),sol2.nrm_v   (:,  1),'-','Color',Colors_trunk{2});
plot(sol1.ToFtrunk+              sol2.tvecbar(:,  1),sol2.nrm_vbar(:,  1),'s','Color',Colors_trunk{2},'MarkerSize',7);
plot(sol1.ToFtrunk+sol2.ToFtrunk+sol2.tvec   (:,end),sol2.nrm_v   (:,end),'-','Color',Colors_branch{ntarg_p_1-2});
plot(sol1.ToFtrunk+sol2.ToFtrunk+sol2.tvecbar(:,end),sol2.nrm_vbar(:,end),'.','Color',Colors_branch{ntarg_p_1-2});

plot(sol1.ToFtrunk+sol2.ToFtrunk+              sol3.tvec   (:,  1),  sol3.nrm_v   (:,  1),  '-','Color',Colors_trunk{3});
plot(sol1.ToFtrunk+sol2.ToFtrunk+              sol3.tvecbar(:,  1),  sol3.nrm_vbar(:,  1),  '^','Color',Colors_trunk{3},'MarkerSize',7);
plot(sol1.ToFtrunk+sol2.ToFtrunk+sol3.ToFtrunk+sol3.tvec   (:,end),  sol3.nrm_v   (:,end),  '-','Color',Colors_branch{ntarg_p_1-3});
plot(sol1.ToFtrunk+sol2.ToFtrunk+sol3.ToFtrunk+sol3.tvecbar(:,end),  sol3.nrm_vbar(:,end),  '.','Color',Colors_branch{ntarg_p_1-3});
plot(sol1.ToFtrunk+sol2.ToFtrunk+sol3.ToFtrunk+sol3.tvec   (:,end-1),sol3.nrm_v   (:,end-1),'-','Color',Colors_branch{ntarg_p_1-4});
plot(sol1.ToFtrunk+sol2.ToFtrunk+sol3.ToFtrunk+sol3.tvecbar(:,end-1),sol3.nrm_vbar(:,end-1),'.','Color',Colors_branch{ntarg_p_1-4});

% title('Speed [m/s]')
ylabel('[m/s]');
xlabel('[s]');
xlim([0,ToFmax])
ylim([0,1.1*prb.vmax]);
ax = gca;
ax.Box = "off";
ax.XLabel.Position = [4.933679040598237,-1.407091365146646,-1];
if saveplot
    exportgraphics(ax,"results/speed_"+interpreter+".pdf",'ContentType','vector');
end

subplot(3,2,3)
plot(linspace(0,ToFmax,100),prb.cost_bound(1)*ones(1,100),'-','Color',Color_bound,'LineWidth',4);
hold on 
plot(              sol1.tvec   (:,  1),sol1.p   (:,  1),'-','Color',Colors_trunk{1});
plot(              sol1.tvecbar(:,  1),sol1.pbar(:,  1),'o','Color',Colors_trunk{1},'MarkerSize',7);
plot(sol1.ToFtrunk+sol1.tvec   (:,end),sol1.p   (:,end),'-','Color',Colors_branch{ntarg_p_1-1});
plot(sol1.ToFtrunk+sol1.tvecbar(:,end),sol1.pbar(:,end),'.','Color',Colors_branch{ntarg_p_1-1});

plot(sol1.ToFtrunk+              sol2.tvec   (:,  1),sol2.p   (:,  1),'-','Color',Colors_trunk{2});
plot(sol1.ToFtrunk+              sol2.tvecbar(:,  1),sol2.pbar(:,  1),'s','Color',Colors_trunk{2},'MarkerSize',7);
plot(sol1.ToFtrunk+sol2.ToFtrunk+sol2.tvec   (:,end),sol2.p   (:,end),'-','Color',Colors_branch{ntarg_p_1-2});
plot(sol1.ToFtrunk+sol2.ToFtrunk+sol2.tvecbar(:,end),sol2.pbar(:,end),'.','Color',Colors_branch{ntarg_p_1-2});

plot(sol1.ToFtrunk+sol2.ToFtrunk+              sol3.tvec   (:,  1),  sol3.p   (:,  1),  '-','Color',Colors_trunk{3});
plot(sol1.ToFtrunk+sol2.ToFtrunk+              sol3.tvecbar(:,  1),  sol3.pbar(:,  1),  '^','Color',Colors_trunk{3},'MarkerSize',7);
plot(sol1.ToFtrunk+sol2.ToFtrunk+sol3.ToFtrunk+sol3.tvec   (:,end),  sol3.p   (:,end),  '-','Color',Colors_branch{ntarg_p_1-3});
plot(sol1.ToFtrunk+sol2.ToFtrunk+sol3.ToFtrunk+sol3.tvecbar(:,end),  sol3.pbar(:,end),  '.','Color',Colors_branch{ntarg_p_1-3});
plot(sol1.ToFtrunk+sol2.ToFtrunk+sol3.ToFtrunk+sol3.tvec   (:,end-1),sol3.p   (:,end-1),'-','Color',Colors_branch{ntarg_p_1-4});
plot(sol1.ToFtrunk+sol2.ToFtrunk+sol3.ToFtrunk+sol3.tvecbar(:,end-1),sol3.pbar(:,end-1),'.','Color',Colors_branch{ntarg_p_1-4});

if interpreter == "latex"
    % title('Trajectory cost [m$^2$/s$^4$]')
    ylabel('[m$^2$/s$^4$]') 
elseif interpreter == "tex"
    % title('Trajectory cost [m^2/s^4]')
    ylabel('[m^2/s^4]')
end
xlabel('[s]');
xlim([0,ToFmax])
ylim([0,1.1*prb.cost_bound(1)]);
ax = gca;
ax.Box = "off";
ax.XLabel.Position = [4.933679040598236,-188.4206573515172,-1];
if saveplot
    exportgraphics(ax,"results/trajcost_"+interpreter+".pdf",'ContentType','vector');
end

subplot(3,2,4)
hold on 
plot(              sol1.tvec   (:,  1),sol1.y   (:,  1),'-','Color',Colors_trunk{1});
plot(              sol1.tvecbar(:,  1),sol1.ybar(:,  1),'o','Color',Colors_trunk{1},'MarkerSize',7);
plot(sol1.ToFtrunk+sol1.tvec   (:,end),sol1.y   (:,end),'-','Color',Colors_branch{ntarg_p_1-1});
plot(sol1.ToFtrunk+sol1.tvecbar(:,end),sol1.ybar(:,end),'.','Color',Colors_branch{ntarg_p_1-1});

plot(sol1.ToFtrunk+              sol2.tvec   (:,  1),sol2.y   (:,  1),'-','Color',Colors_trunk{2});
plot(sol1.ToFtrunk+              sol2.tvecbar(:,  1),sol2.ybar(:,  1),'s','Color',Colors_trunk{2},'MarkerSize',7);
plot(sol1.ToFtrunk+sol2.ToFtrunk+sol2.tvec   (:,end),sol2.y   (:,end),'-','Color',Colors_branch{ntarg_p_1-2});
plot(sol1.ToFtrunk+sol2.ToFtrunk+sol2.tvecbar(:,end),sol2.ybar(:,end),'.','Color',Colors_branch{ntarg_p_1-2});

plot(sol1.ToFtrunk+sol2.ToFtrunk+              sol3.tvec   (:,  1),  sol3.y   (:,  1),  '-','Color',Colors_trunk{3});
plot(sol1.ToFtrunk+sol2.ToFtrunk+              sol3.tvecbar(:,  1),  sol3.ybar(:,  1),  '^','Color',Colors_trunk{3},'MarkerSize',7);
plot(sol1.ToFtrunk+sol2.ToFtrunk+sol3.ToFtrunk+sol3.tvec   (:,end),  sol3.y   (:,end),  '-','Color',Colors_branch{ntarg_p_1-3});
plot(sol1.ToFtrunk+sol2.ToFtrunk+sol3.ToFtrunk+sol3.tvecbar(:,end),  sol3.ybar(:,end),  '.','Color',Colors_branch{ntarg_p_1-3});
plot(sol1.ToFtrunk+sol2.ToFtrunk+sol3.ToFtrunk+sol3.tvec   (:,end-1),sol3.y   (:,end-1),'-','Color',Colors_branch{ntarg_p_1-4});
plot(sol1.ToFtrunk+sol2.ToFtrunk+sol3.ToFtrunk+sol3.tvecbar(:,end-1),sol3.ybar(:,end-1),'.','Color',Colors_branch{ntarg_p_1-4});

% title('Constraint violation')
xlabel('[s]');
xlim([0,ToFmax]);
ax = gca;
ax.Box = "off";
ax.XLabel.Position = [4.933679040598238,-0.000003398020972,-1];
if saveplot
    exportgraphics(ax,"results/cnstrviol_"+interpreter+".pdf",'ContentType','vector');
end

subplot(3,2,5)
plot(linspace(0,ToFmax,100),prb.Tmax*ones(1,100),'-','Color',Color_bound,'LineWidth',4);
hold on 
plot(linspace(0,ToFmax,100),prb.Tmin*ones(1,100),'-','Color',Color_bound,'LineWidth',4);
plot(              sol1.tvec   (:,  1),sol1.nrm_T   (:,  1),'-','Color',Colors_trunk{1});
plot(              sol1.tvecbar(:,  1),sol1.nrm_Tbar(:,  1),'o','Color',Colors_trunk{1},'MarkerSize',7);
plot(sol1.ToFtrunk+sol1.tvec   (:,end),sol1.nrm_T   (:,end),'-','Color',Colors_branch{ntarg_p_1-1});
plot(sol1.ToFtrunk+sol1.tvecbar(:,end),sol1.nrm_Tbar(:,end),'.','Color',Colors_branch{ntarg_p_1-1});

plot(sol1.ToFtrunk+              sol2.tvec   (:,  1),sol2.nrm_T   (:,  1),'-','Color',Colors_trunk{2});
plot(sol1.ToFtrunk+              sol2.tvecbar(:,  1),sol2.nrm_Tbar(:,  1),'s','Color',Colors_trunk{2},'MarkerSize',7);
plot(sol1.ToFtrunk+sol2.ToFtrunk+sol2.tvec   (:,end),sol2.nrm_T   (:,end),'-','Color',Colors_branch{ntarg_p_1-2});
plot(sol1.ToFtrunk+sol2.ToFtrunk+sol2.tvecbar(:,end),sol2.nrm_Tbar(:,end),'.','Color',Colors_branch{ntarg_p_1-2});

plot(sol1.ToFtrunk+sol2.ToFtrunk+              sol3.tvec   (:,  1),  sol3.nrm_T   (:,  1),  '-','Color',Colors_trunk{3});
plot(sol1.ToFtrunk+sol2.ToFtrunk+              sol3.tvecbar(:,  1),  sol3.nrm_Tbar(:,  1),  '^','Color',Colors_trunk{3},'MarkerSize',7);
plot(sol1.ToFtrunk+sol2.ToFtrunk+sol3.ToFtrunk+sol3.tvec   (:,end),  sol3.nrm_T   (:,end),  '-','Color',Colors_branch{ntarg_p_1-3});
plot(sol1.ToFtrunk+sol2.ToFtrunk+sol3.ToFtrunk+sol3.tvecbar(:,end),  sol3.nrm_Tbar(:,end),  '.','Color',Colors_branch{ntarg_p_1-3});
plot(sol1.ToFtrunk+sol2.ToFtrunk+sol3.ToFtrunk+sol3.tvec   (:,end-1),sol3.nrm_T   (:,end-1),'-','Color',Colors_branch{ntarg_p_1-4});
plot(sol1.ToFtrunk+sol2.ToFtrunk+sol3.ToFtrunk+sol3.tvecbar(:,end-1),sol3.nrm_Tbar(:,end-1),'.','Color',Colors_branch{ntarg_p_1-4});

if interpreter == "latex"
    % title('Thrust magnitude [m/s$^2$]')
    ylabel('[m/s$^2$]');
elseif interpreter == "tex"
    % title('Thrust magnitude [m/s^2]')
    ylabel('[m/s^2]');
end
xlabel('[s]');
xlim([0,ToFmax])
ylim([0,1.1*prb.Tmax]);
ax = gca;
ax.Box = "off";
ax.XLabel.Position = [4.933679040598236,-3.566516881046889,-0.999999999999986];
if saveplot
    exportgraphics(ax,"results/thrust_"+interpreter+".pdf",'ContentType','vector');
end

savefig(f2,"results/data");