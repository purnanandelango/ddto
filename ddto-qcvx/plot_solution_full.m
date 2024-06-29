clc
close all
clearvars

set_params

load solution

ntarg = length(target_order);

interpreter = "tex";        plt.setfig("tex");
% interpreter = "latex";      plt.setfig("latex");
saveplot = true;

Colors_branch = {[1,0.5,0.5],[0,0.7,0.3],[0.5,0.5,1],[1,0.65,0.3]};
Colors_trunk = {[0.1,0.1,0.2],[0.3,0.3,0.4],[0.5,0.5,0.6]};
Color_bound = [161,240,242]/255;

f1 = figure('Position',[100,100,700,700]);
hold on

% Plot trunk and branch position trace
idx = find(prob_par.tag==target_order(1));
[~,N] = size(soldd.X{idx});

plot3(soldd.X{1}(1,1:soldd.dd_idx+1),soldd.X{1}(2,1:soldd.dd_idx+1),soldd.X{1}(3,1:soldd.dd_idx+1),'o-','Color',Colors_trunk{1},'MarkerSize',7);
plot3(soldd.X{idx}(1,soldd.dd_idx+1:N),soldd.X{idx}(2,soldd.dd_idx+1:N),soldd.X{idx}(3,soldd.dd_idx+1:N),'.-','Color',Colors_branch{ntarg});

idx = find(prob_par2.tag==target_order(2));
[~,N] = size(soldd2.X{idx});

plot3(soldd2.X{1}(1,1:soldd2.dd_idx+1),soldd2.X{1}(2,1:soldd2.dd_idx+1),soldd2.X{1}(3,1:soldd2.dd_idx+1),'s-','Color',Colors_trunk{2},'MarkerSize',7);
plot3(soldd2.X{idx}(1,soldd2.dd_idx+1:N),soldd2.X{idx}(2,soldd2.dd_idx+1:N),soldd2.X{idx}(3,soldd2.dd_idx+1:N),'.-','Color',Colors_branch{ntarg-1});

idx = find(prob_par3.tag==target_order(3));
[~,N] = size(soldd3.X{idx});

plot3(soldd3.X{1}(1,1:soldd3.dd_idx+1),soldd3.X{1}(2,1:soldd3.dd_idx+1),soldd3.X{1}(3,1:soldd3.dd_idx+1),'^-','Color',Colors_trunk{3},'MarkerSize',7);
plot3(soldd3.X{idx}(1,soldd3.dd_idx+1:N),soldd3.X{idx}(2,soldd3.dd_idx+1:N),soldd3.X{idx}(3,soldd3.dd_idx+1:N),'.-','Color',Colors_branch{ntarg-2});

idx = find(prob_par3.tag==target_order(4));
[~,N] = size(soldd3.X{idx});

plot3(soldd3.X{idx}(1,soldd3.dd_idx+1:N),soldd3.X{idx}(2,soldd3.dd_idx+1:N),soldd3.X{idx}(3,soldd3.dd_idx+1:N),'.-','Color',Colors_branch{ntarg-3});

grid on

ax = gca;
ax.Box = "off";
ax.View = [68,10];
ax.ZLim = [-3,33];
ax.XLim = [-1,41];
ax.YLim = [-7,41];
ax.XLabel.String = '[m]';
ax.XLabel.Position = [19.78481709803691,-9.502394273346567,-7.184835973095858];
ax.YLabel.String = '[m]';
ax.YLabel.Position = [42.114198364688946,17.489138153402337,-7.819637146078605];
ax.ZLabel.String = '[m]';

ax.DataAspectRatio = [1,1,1];
ax.PlotBoxAspectRatio = [1,1,1];

if interpreter == "tex"
    annotation("textbox",'Position',[0.679999999999999,0.312571428571427,0.088571428571429,0.076428571428571],'EdgeColor','none','String','{\itz}^4');
    annotation("textbox",'Position',[0.68,0.381142857142857,0.088571428571429,0.076428571428571],'EdgeColor','none','String','{\itz}^3');
    annotation("textbox",'Position',[0.47,0.271142857142857,0.088571428571429,0.076428571428571],'EdgeColor','none','String','{\itz}^2');
    annotation("textbox",'Position',[0.32,0.259714285714286,0.088571428571429,0.076428571428571],'EdgeColor','none','String','{\itz}^1');
    annotation("textbox",'Position',[0.201428571428571,0.728285714285714,0.088571428571429,0.076428571428571],'EdgeColor','none','String','{\itz}^0');
elseif interpreter == "latex"
    annotation("textbox",'Position',[0.678571428571428,0.296857142857141,0.088571428571429,0.076428571428571],'EdgeColor','none','String','$z^4$');
    annotation("textbox",'Position',[0.681428571428571,0.366857142857143,0.088571428571429,0.076428571428571],'EdgeColor','none','String','$z^3$');
    annotation("textbox",'Position',[0.468571428571429,0.255428571428572,0.088571428571429,0.076428571428571],'EdgeColor','none','String','$z^2$');
    annotation("textbox",'Position',[0.318571428571429,0.245428571428572,0.088571428571429,0.076428571428571],'EdgeColor','none','String','$z^1$');
    annotation("textbox",'Position',[0.201428571428571,0.714,0.088571428571429,0.076428571428571],'EdgeColor','none','String','$z^0$');
end

if saveplot 
    exportgraphics(ax,"results/position_"+interpreter+".pdf",'ContentType','vector');
    savefig(f1,'results/position.fig');    
end

Nmax = max(prob_par.N);

f2 = figure('Position',[1500,101,1400,700],'Units','normalized');

% Thrust pointing angle
subplot(3,2,1)
plot(linspace(0,Nmax,100),atand(alf_tp)*ones(1,100),'-','Color',Color_bound,'LineWidth',4);
hold on

idx = find(prob_par.tag==target_order(1));
[~,N] = size(soldd.X{idx});

point_angle = acosd(soldd.U{1}(3,1:soldd.dd_idx)./vecnorm(soldd.U{1}(:,1:soldd.dd_idx)));
point_angle(end+1) = point_angle(end);
stairs(1:soldd.dd_idx+1,point_angle,'o-','Color',Colors_trunk{1},'MarkerSize',7);
point_angle = acosd(soldd.U{idx}(3,soldd.dd_idx+1:N-1)./vecnorm(soldd.U{idx}(:,soldd.dd_idx+1:N-1)));
point_angle(end+1) = point_angle(end);
stairs(soldd.dd_idx+1:N,point_angle,'.-','Color',Colors_branch{ntarg},'MarkerSize',20);

idx = find(prob_par2.tag==target_order(2));
[~,N] = size(soldd2.X{idx});

if soldd2.dd_idx > 0
point_angle = acosd(soldd2.U{1}(3,1:soldd2.dd_idx)./vecnorm(soldd2.U{1}(:,1:soldd2.dd_idx)));
point_angle(end+1) = point_angle(end);
stairs(soldd.dd_idx+(1:soldd2.dd_idx+1),point_angle,'s-','Color',Colors_trunk{2},'MarkerSize',7);
end
point_angle = acosd(soldd2.U{idx}(3,soldd2.dd_idx+1:N-1)./vecnorm(soldd2.U{idx}(:,soldd2.dd_idx+1:N-1)));
point_angle(end+1) = point_angle(end);
stairs(soldd.dd_idx+(soldd2.dd_idx+1:N),point_angle,'.-','Color',Colors_branch{ntarg-1},'MarkerSize',20);

idx = find(prob_par3.tag==target_order(3));
[~,N] = size(soldd3.X{idx});

point_angle = acosd(soldd3.U{1}(3,1:soldd3.dd_idx)./vecnorm(soldd3.U{1}(:,1:soldd3.dd_idx)));
point_angle(end+1) = point_angle(end);
stairs(soldd.dd_idx+soldd2.dd_idx+(1:soldd3.dd_idx+1),point_angle,'^-','Color',Colors_trunk{3},'MarkerSize',7);
point_angle = acosd(soldd3.U{idx}(3,soldd3.dd_idx+1:N-1)./vecnorm(soldd3.U{idx}(:,soldd3.dd_idx+1:N-1)));
point_angle(end+1) = point_angle(end);
stairs(soldd.dd_idx+soldd2.dd_idx+(soldd3.dd_idx+1:N),point_angle,'.-','Color',Colors_branch{ntarg-2},'MarkerSize',20);

idx = find(prob_par3.tag==target_order(4));
[~,N] = size(soldd3.X{idx});
point_angle = acosd(soldd3.U{idx}(3,soldd3.dd_idx+1:N-1)./vecnorm(soldd3.U{idx}(:,soldd3.dd_idx+1:N-1)));
point_angle(end+1) = point_angle(end);
stairs(soldd.dd_idx+soldd2.dd_idx+(soldd3.dd_idx+1:N),point_angle,'.-','Color',Colors_branch{ntarg-3},'MarkerSize',20);

ylabel('[deg]');
xlabel('[s]');
xlim([1,Nmax]);
xticks([1,4,8,12,16,20]);
ylim([0,1.1*atand(alf_tp)]);
ax = gca;
ax.Box = "off";
ax.XLabel.Position = [10.418899134612507,-9.469637930063676,-1];
if saveplot
    exportgraphics(ax,"results/pointangle_"+interpreter+".pdf");
end

% Cost bound
cost_bound = 2*max(solopt.cost);

subplot(3,2,2)
plot(linspace(0,Nmax,100),cost_bound*ones(1,100),'-','Color',Color_bound,'LineWidth',4);
hold on

idx = find(prob_par.tag==target_order(1));
[~,N] = size(soldd.X{idx});

cumcost1 = cumsum(sum(soldd.U{1}(:,1:soldd.dd_idx) .^ 2,1));
cumcost1(end+1) = cumcost1(end);
stairs(1:soldd.dd_idx+1,cumcost1,'o-','Color',Colors_trunk{1},'MarkerSize',7);
cumcost = cumsum(sum(soldd.U{idx}(:,soldd.dd_idx+1:N-1) .^ 2,1));
cumcost(end+1) = cumcost(end);
stairs(soldd.dd_idx+1:N,cumcost1(end)+cumcost,'.-','Color',Colors_branch{ntarg},'MarkerSize',20);

idx = find(prob_par2.tag==target_order(2));
[~,N] = size(soldd2.X{idx});

if soldd2.dd_idx > 0
cumcost2 = cumsum(sum(soldd2.U{1}(:,1:soldd2.dd_idx) .^ 2,1));
cumcost2(end+1) = cumcost2(end);
else
cumcost2 = 0;
end
stairs(soldd.dd_idx+(1:soldd2.dd_idx+1),cumcost1(end)+cumcost2,'s-','Color',Colors_trunk{2},'MarkerSize',7);
cumcost = cumsum(sum(soldd2.U{idx}(:,soldd2.dd_idx+1:N-1) .^ 2,1));
cumcost(end+1) = cumcost(end);
stairs(soldd.dd_idx+(soldd2.dd_idx+1:N),cumcost1(end)+cumcost2(end)+cumcost,'.-','Color',Colors_branch{ntarg-1},'MarkerSize',20);

idx = find(prob_par3.tag==target_order(3));
[~,N] = size(soldd3.X{idx});

cumcost3 = cumsum(sum(soldd3.U{1}(:,1:soldd3.dd_idx) .^ 2,1));
cumcost3(end+1) = cumcost3(end);
stairs(soldd.dd_idx+soldd2.dd_idx+(1:soldd3.dd_idx+1),cumcost1(end)+cumcost2(end)+cumcost3,'^-','Color',Colors_trunk{3},'MarkerSize',7);
cumcost = cumsum(sum(soldd3.U{idx}(:,soldd3.dd_idx+1:N-1) .^ 2,1));
cumcost(end+1) = cumcost(end);
stairs(soldd.dd_idx+soldd2.dd_idx+(soldd3.dd_idx+1:N),cumcost1(end)+cumcost2(end)+cumcost3(end)+cumcost,'.-','Color',Colors_branch{ntarg-2},'MarkerSize',20);

idx = find(prob_par3.tag==target_order(4));
[~,N] = size(soldd3.X{idx});
cumcost = cumsum(sum(soldd3.U{idx}(:,soldd3.dd_idx+1:N-1) .^ 2,1));
cumcost(end+1) = cumcost(end);
stairs(soldd.dd_idx+soldd2.dd_idx+(soldd3.dd_idx+1:N),cumcost1(end)+cumcost2(end)+cumcost3(end)+cumcost,'.-','Color',Colors_branch{ntarg-3},'MarkerSize',20);

if interpreter == "tex"
    ylabel('[m^2/s^4]');
elseif interpreter == "latex"
    ylabel('[m$^2$/s$^4$]');
end
xlabel('[s]');
xlim([1,Nmax]);
xticks([1,4,8,12,16,20]);
ylim([0,1.1*cost_bound]);
ax = gca;
ax.Box = "off";
ax.XLabel.Position = [10.500009059906006,-581.1751560117887,-1];
if saveplot
    exportgraphics(ax,"results/trajcost_"+interpreter+".pdf");
end

% Thrust magnitude
subplot(3,2,3)
hold on
plot(linspace(0,Nmax,100),u_max*ones(1,100),'-','Color',Color_bound,'LineWidth',4);
plot(linspace(0,Nmax,100),u_min*ones(1,100),'-','Color',Color_bound,'LineWidth',4);

idx = find(prob_par.tag==target_order(1));
[~,N] = size(soldd.X{idx});

nrm_T = vecnorm(soldd.U{1}(:,1:soldd.dd_idx));
nrm_T(end+1) = nrm_T(end);
stairs(1:soldd.dd_idx+1,nrm_T,'o-','Color',Colors_trunk{1},'MarkerSize',7);
nrm_T = vecnorm(soldd.U{idx}(:,soldd.dd_idx+1:N-1));
nrm_T(end+1) = nrm_T(end);
stairs(soldd.dd_idx+1:N,nrm_T,'.-','Color',Colors_branch{ntarg},'MarkerSize',20);

idx = find(prob_par2.tag==target_order(2));
[~,N] = size(soldd2.X{idx});

if soldd2.dd_idx > 0
nrm_T = vecnorm(soldd2.U{1}(:,1:soldd2.dd_idx));
nrm_T(end+1) = nrm_T(end);
stairs(soldd.dd_idx+(1:soldd2.dd_idx+1),nrm_T,'s-','Color',Colors_trunk{2},'MarkerSize',7);
end
nrm_T = vecnorm(soldd2.U{idx}(:,soldd2.dd_idx+1:N-1));
nrm_T(end+1) = nrm_T(end);
stairs(soldd.dd_idx+(soldd2.dd_idx+1:N),nrm_T,'.-','Color',Colors_branch{ntarg-1},'MarkerSize',20);

idx = find(prob_par3.tag==target_order(3));
[~,N] = size(soldd3.X{idx});

nrm_T = vecnorm(soldd3.U{1}(:,1:soldd3.dd_idx));
nrm_T(end+1) = nrm_T(end);
stairs(soldd.dd_idx+soldd2.dd_idx+(1:soldd3.dd_idx+1),nrm_T,'^-','Color',Colors_trunk{3},'MarkerSize',7);
nrm_T = vecnorm(soldd3.U{idx}(:,soldd3.dd_idx+1:N-1));
nrm_T(end+1) = nrm_T(end);
stairs(soldd.dd_idx+soldd2.dd_idx+(soldd3.dd_idx+1:N),nrm_T,'.-','Color',Colors_branch{ntarg-2},'MarkerSize',20);

idx = find(prob_par3.tag==target_order(4));
[~,N] = size(soldd3.X{idx});
nrm_T = vecnorm(soldd3.U{idx}(:,soldd3.dd_idx+1:N-1));
nrm_T(end+1) = nrm_T(end);
stairs(soldd.dd_idx+soldd2.dd_idx+(soldd3.dd_idx+1:N),nrm_T,'.-','Color',Colors_branch{ntarg-3},'MarkerSize',20);

if interpreter == "tex"
    ylabel('[m/s^2]');
elseif interpreter == "latex"
    ylabel('[m/s$^2$]');
end
xlabel('[s]');
xlim([1,Nmax]);
xticks([1,4,8,12,16,20]);
ylim([0,1.1*u_max]);
ax = gca;
ax.Box = "off";
ax.XLabel.Position = [10.418899134612511,-3.037744242356034,-0.999999999999986];
if saveplot
    exportgraphics(ax,"results/thrust_"+interpreter+".pdf");
end

savefig(f2,"results/data");