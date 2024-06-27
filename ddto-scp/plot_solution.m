clearvars
close all

load recent_solution.mat

ntarg_p_1 = prb.ntarg+1;
ToFtrunk = tvecbar(end,1);
ToFmax = max(ToFtrunk+tvec(end,2:end));

figure

saveplot = false;

Colors = {[1,0.5,0.5],[0.5,1,0.5],[0.5,0.5,1],[1,1,0.5],[1,0.5,1],[0.5,1,1]};

% Plot trunk position
plot3(r(1,:,1),r(2,:,1),r(3,:,1),'-k');
hold on 
plot3(rbar(1,:,1),rbar(2,:,1),rbar(3,:,1),'.k');
for j = 2:ntarg_p_1
    plot3(r(1,:,j),r(2,:,j),r(3,:,j),'-','Color',Colors{j});
    hold on 
    plot3(rbar(1,:,j),rbar(2,:,j),rbar(3,:,j),'.','Color',Colors{j});
end

grid on

ax = gca;
ax.DataAspectRatio = [1,1,1];
ax.PlotBoxAspectRatio = [1,1,1];    

[X1,Y1,Z1] = ellipsoid(0,0,0,1/0.3,1/0.1,1/0.3,100);
X1 = X1 + prb.qobs(1,1);
Y1 = Y1 + prb.qobs(2,1);
Z1 = Z1 + prb.qobs(3,1);
ellip1 = surf(X1,Y1,Z1,'EdgeColor',[1,0,0],'EdgeAlpha',0.5,'FaceColor',[1,0,0],'FaceAlpha',0.2);
rotate(ellip1,[0,0,1],0);

[X2,Y2,Z2] = ellipsoid(0,0,0,1/0.1,1/0.3,1/0.3,100);
X2 = X2 + prb.qobs(1,2);
Y2 = Y2 + prb.qobs(2,2);
Z2 = Z2 + prb.qobs(3,2);
ellip2 = surf(X2,Y2,Z2,'EdgeColor',[0,0,1],'EdgeAlpha',0.5,'FaceColor',[0,0,1],'FaceAlpha',0.2);
rotate(ellip2,[0,0,1],0);    

view(-6,1);
title('Position [m]');

if saveplot
    exportgraphics(ax,'results/position.pdf','BackgroundColor','none');
end

figure('Position',[100,100,1600,1000])

point_angle(prb.Kfine,ntarg_p_1) = 0;
nrm_T(prb.Kfine,ntarg_p_1) = 0;
nrm_v(prb.Kfine,ntarg_p_1) = 0;
point_anglebar(prb.K,ntarg_p_1) = 0;
nrm_Tbar(prb.K,ntarg_p_1) = 0;
nrm_vbar(prb.K,ntarg_p_1) = 0;
for k = 1:prb.Kfine
    for j = 1:ntarg_p_1
        point_angle(k,j) = acosd((prb.ehat'*T(:,k,j))/norm(T(:,k,j)));
        nrm_T(k,j) = norm(T(:,k,j));
        nrm_v(k,j) = norm(v(:,k,j));
        if k <= prb.K
            point_anglebar(k,j) = acosd((prb.ehat'*Tbar(:,k,j))/norm(Tbar(:,k,j)));
            nrm_Tbar(k,j) = norm(Tbar(:,k,j));
            nrm_vbar(k,j) = norm(vbar(:,k,j));
        end
    end
end

subplot(3,2,1)
plot(linspace(0,ToFmax,100),prb.deltamax*ones(1,100),'-r','LineWidth',4);
hold on 
plot(tvec(:,1),point_angle(:,1),'-k');
plot(tvecbar(:,1),point_anglebar(:,1),'.k');
for j = 2:ntarg_p_1
    plot(ToFtrunk+tvec(:,j),point_angle(:,j),'-','Color',Colors{j});
    plot(ToFtrunk+tvecbar(:,j),point_anglebar(:,j),'.','Color',Colors{j});    
end
title('Thrust pointing angle [deg]')
xlabel('$t$ [s]');
xlim([0,ToFmax])
ylim([0,1.1*prb.deltamax]);
if saveplot
    ax = gca;
    exportgraphics(ax,'results/pointangle.pdf','BackgroundColor','none');
end

subplot(3,2,2)
plot(linspace(0,ToFmax,100),prb.vmax*ones(1,100),'-r','LineWidth',4);
hold on 
plot(tvec(:,1),nrm_v(:,1),'-k');
plot(tvecbar(:,1),nrm_vbar(:,1),'.k');
for j = 2:ntarg_p_1
    plot(ToFtrunk+tvec(:,j),nrm_v(:,j),'-','Color',Colors{j});
    plot(ToFtrunk+tvecbar(:,j),nrm_vbar(:,j),'.','Color',Colors{j});    
end
title('Speed [m/s]')
xlabel('$t$ [s]');
xlim([0,ToFmax])
ylim([0,1.1*prb.vmax]);
if saveplot
    ax = gca;
    exportgraphics(ax,'results/speed.pdf','BackgroundColor','none');
end

subplot(3,2,3)
plot(tvec(:,1),p(:,1),'-k');
hold on
plot(tvecbar(:,1),pbar(:,1),'.k');
for j = 2:ntarg_p_1
    plot(linspace(0,ToFmax,100),prb.cost_bound(j-1)*ones(1,100),'.','MarkerSize',4,'Color',Colors{j});
    plot(ToFtrunk+tvec(:,j),p(:,j),'-','Color',Colors{j});
    plot(ToFtrunk+tvecbar(:,j),pbar(:,j),'.','Color',Colors{j});    
end
title('Trajectory cost [m$^2$/s$^{4}$]')
xlabel('$t$ [s]');
xlim([0,ToFmax])
ylim([0,1.1*prb.pmax]);
if saveplot
    ax = gca;
    exportgraphics(ax,'results/trajcost.pdf','BackgroundColor','none');
end

subplot(3,2,4)
plot(tvec(:,1),y(:,1),'-k');
hold on
plot(tvecbar(:,1),ybar(:,1),'.k');
for j = 2:ntarg_p_1
    plot(ToFtrunk+tvec(:,j),y(:,j),'-','Color',Colors{j});
    plot(ToFtrunk+tvecbar(:,j),ybar(:,j),'.','Color',Colors{j});    
end
title('Constraint violation')
xlabel('$t$ [s]');
xlim([0,ToFmax])
if saveplot
    ax = gca;
    exportgraphics(ax,'results/cnstrviol.pdf','BackgroundColor','none');
end

subplot(3,2,5)
plot(linspace(0,ToFmax,100),prb.Tmax*ones(1,100),'-r','LineWidth',4);
hold on 
plot(tvec(:,1),nrm_T(:,1),'-k');
plot(tvecbar(:,1),nrm_Tbar(:,1),'.k');
for j = 2:ntarg_p_1
    plot(ToFtrunk+tvec(:,j),nrm_T(:,j),'-','Color',Colors{j});
    plot(ToFtrunk+tvecbar(:,j),nrm_Tbar(:,j),'.','Color',Colors{j});    
end
title('Thrust magnitude [m/s$^2$]')
xlabel('$t$ [s]');
xlim([0,ToFmax])
ylim([0,1.1*prb.Tmax]);
if saveplot
    ax = gca;
    exportgraphics(ax,'results/thrust.pdf','BackgroundColor','none');
end

% subplot(3,2,6)
% plot(prb.tau,tvecbar(:,1),'.k');
% hold on
% plot(tau,tvec(:,1),'-k');
% for j=2:ntarg_p_1
%     plot(prb.tau,ToFtrunk+tvecbar(:,j),'.','Color',Colors{j});
%     plot(tau,ToFtrunk+tvec(:,j),'-','Color',Colors{j});
% end
% xlabel('$\tau$');
% ylim([0,1.1*ToFmax]);
% title('Time [s]')
% if saveplot
%     ax = gca;
%     exportgraphics(ax,'results/time.pdf','BackgroundColor','none');
% end

subplot(3,2,6)
plot(prb.tau,prb.smin*ones(1,prb.K),'-r','LineWidth',4);
hold on
plot(prb.tau,prb.smax*ones(1,prb.K),'-r','LineWidth',4);
plot(prb.tau,sbar(:,1),'.k');
p2 = plot(tau,s(:,1),'-k');
for j=2:ntarg_p_1
    plot(prb.tau,sbar(:,j),'.','Color',Colors{j});
    plot(tau,s(:,j),'-','Color',Colors{j});
end
xlabel('$\tau$');
ylim([0,1.1*prb.smax]);
title('Dilation factor')
if saveplot
    ax = gca;
    exportgraphics(ax,'results/dilation.pdf','BackgroundColor','none');
end