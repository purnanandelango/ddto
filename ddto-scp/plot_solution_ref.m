clearvars
close all

load recent_solution.mat

figure

saveplot = false;

suffix = '3D';

plot3(r(1,:),r(2,:),r(3,:),'-k');
hold on 
plot3(xbar(1,:),xbar(2,:),xbar(3,:),'.k');
grid on

ax = gca;
ax.DataAspectRatio = [1,1,1];
ax.PlotBoxAspectRatio = [1,1,1];    

[X1,Y1,Z1] = ellipsoid(0,0,0,1/0.3,1/0.1,1/0.3);
X1 = X1 + prb.qobs(1,1);
Y1 = Y1 + prb.qobs(2,1);
Z1 = Z1 + prb.qobs(3,1);
ellip1 = surf(X1,Y1,Z1,'EdgeColor',[1,0,0],'EdgeAlpha',0.5,'FaceColor',[1,0,0],'FaceAlpha',0.2);
rotate(ellip1,[0,0,1],0);

[X2,Y2,Z2] = ellipsoid(0,0,0,1/0.1,1/0.3,1/0.3);
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

figure
point_angle(prb.Kfine) = 0;
nrm_T(prb.Kfine) = 0;
nrm_v(prb.Kfine) = 0;
for j = 1:prb.Kfine
    point_angle(j) = acosd((prb.ehat'*u(1:prb.n,j))/norm(u(1:prb.n,j)));
    nrm_T(j) = norm(u(1:prb.n,j));
    nrm_v(j) = norm(v(1:prb.n,j));
end
point_anglebar(prb.K) = 0;
nrm_Tbar(prb.K) = 0;
nrm_vbar(prb.K) = 0;
for j = 1:prb.K
    point_anglebar(j) = acosd((prb.ehat'*ubar(1:prb.n,j))/norm(ubar(1:prb.n,j)));
    nrm_Tbar(j) = norm(ubar(1:prb.n,j));
    nrm_vbar(j) = norm(xbar(prb.n+1:2*prb.n,j));
end

subplot(2,2,1)
plot(tvecbar,prb.deltamax*ones(1,prb.K),'-r');
hold on 
plot(tvec,point_angle,'-b');
plot(tvecbar,point_anglebar,'.b');
title('Pointing angle [deg]')
xlabel('$t$ [s]');
xlim([0,tvec(end)])
ylim([0,1.1*prb.deltamax]);
if saveplot
    ax = gca;
    exportgraphics(ax,'results/pointangle.pdf','BackgroundColor','none');
end

subplot(2,2,2)
plot(tvecbar,prb.vmax*ones(1,prb.K),'-r');
hold on 
plot(tvec,nrm_v,'-b');
plot(tvecbar,nrm_vbar,'.b');
title('Speed [m s$^{-1}$]')
xlabel('$t$ [s]');
xlim([0,tvec(end)])
ylim([0,1.1*prb.vmax]);
if saveplot
    ax = gca;
    exportgraphics(ax,'results/speed.pdf','BackgroundColor','none');
end

subplot(2,2,3)
plot(tvecbar,prb.Tmin*ones(1,prb.K),'-r');
hold on 
plot(tvecbar,prb.Tmax*ones(1,prb.K),'-r');
plot(tvec,nrm_T,'-b');
plot(tvecbar,nrm_Tbar,'.b');
title('Thrust [m s$^{-2}$]');
xlabel('$t$ [s]');
xlim([0,tvec(end)])
ylim([0,1.1*prb.Tmax])
if saveplot
    ax = gca;
    exportgraphics(ax,'results/thrust.pdf','BackgroundColor','none');
end

subplot(2,2,4)
plot(prb.tau,prb.smin*ones(1,prb.K),'-r');
hold on
plot(prb.tau,prb.smax*ones(1,prb.K),'-r');
plot(prb.tau,tvecbar,'.k');
p1 = plot(tau,tvec,'-k');
plot(prb.tau,ubar(prb.n+1,:),'.g');
p2 = plot(tau,u(prb.n+1,:),'-g');
legend([p1,p2],{'$t(\tau)$','$s(\tau)$'})
xlabel('$\tau$');
ylim([0,1.1*prb.smax]);
title('Time \& Dilation')
if saveplot
    ax = gca;
    exportgraphics(ax,'results/time_dilation.pdf','BackgroundColor','none');
end