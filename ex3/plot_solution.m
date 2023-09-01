clearvars
close all

load recent_solution
cval = {[0.6,0.4,0],[0.6,0,0.4],[0,0.4,0.6],[0.6,0.4,0.6]};

figure
subplot(2,2,1)
plot(r(1,:,1),r(2,:,1),'-k');
hold on 
plot(rbar(1,:,1),rbar(2,:,1),'+k');    
for j = 2:prb.ntarg+1        
    plot(r(1,:,j),r(2,:,j),'-','Color',cval{j-1});
    plot(rbar(1,:,j),rbar(2,:,j),'o','Color',cval{j-1});
end
title('Position');

th = linspace(0,2*pi);
for j = 1:prb.nobs
    pobs = prb.robs(:,j) + prb.aobs(j)*[cos(th);sin(th)];
    plot(pobs(1,:),pobs(2,:),'-k','LineWidth',0.5);
end
xlim([-5,6]);
ylim([-3,8]);

ax = gca;
ax.DataAspectRatio = [1,1,1];
ax.PlotBoxAspectRatio = [1,1,1];


subplot(2,2,2)
plot(tvec(:,1),nrmv(:,1),'-k')
hold on
plot(tvecbar(:,1),nrmvbar(:,1),'ok');
plot(tvecbar(:,1),prb.vmax*ones(1,prb.K),'-r','LineWidth',2);
for j = 2:prb.ntarg+1    
    plot(tvecbar(end,1)+tvecbar(:,j),prb.vmax*ones(1,prb.K),'-r','LineWidth',2);    
    plot(tvec(end,1)+tvec(:,j),nrmv(:,j),'-','Color',cval{j-1});
    plot(tvecbar(end,1)+tvecbar(:,j),nrmvbar(:,j),'o','Color',cval{j-1});
end
xlim([0,tvecbar(end,1)+max(tvecbar(end,2:end))]);
title('Speed')
xlabel('$t$');

subplot(2,2,3)
plot(tvec(:,1),nrmu(:,1),'-k');
hold on
plot(tvecbar(:,1),nrmubar(:,1),'ok');    
plot(tvecbar(:,1),prb.umax*ones(1,prb.K),'-r','LineWidth',2);
plot(tvecbar(:,1),prb.umin*ones(1,prb.K),'-r','LineWidth',2);    
for j = 2:prb.ntarg+1
    plot(tvecbar(end,1)+tvecbar(:,j),prb.umax*ones(1,prb.K),'-r','LineWidth',2);
    plot(tvecbar(end,1)+tvecbar(:,j),prb.umin*ones(1,prb.K),'-r','LineWidth',2);      
    plot(tvec(end,1)+tvec(:,j),nrmu(:,j),'-','Color',cval{j-1});
    plot(tvecbar(end,1)+tvecbar(:,j),nrmubar(:,j),'o','Color',cval{j-1});
end
xlim([0,tvecbar(end,1)+max(tvecbar(end,2:end))]);
title('Thrust magnitude');
xlabel('$t$');

subplot(2,2,4)
plt_t = plot(prb.tau,tvecbar(:,1),'kd');
hold on
plot(tau,tvec(:,1),'-k');    
plt_dt = plot(prb.tau(1:end-1),diff(tvecbar(:,1)),'-ks');    
plt_s = plot(prb.tau,sbar(:,1),'k+');
plot(tau,s(:,1),'-k');
for j = 2:prb.ntarg+1
    plt_t = plot(prb.tau,tvecbar(:,j),'d','Color',cval{j-1});
    hold on
    plot(tau,tvec(:,j),'-','Color',cval{j-1});    
    plt_dt = plot(prb.tau(1:end-1),diff(tvecbar(:,j)),'-s','Color',cval{j-1});    
    plt_s = plot(prb.tau,sbar(:,j),'+','Color',cval{j-1});
    plot(tau,s(:,j),'-','Color',cval{j-1});    
end
legend([plt_t,plt_dt,plt_s],{'$t$','$\Delta t$','$s$'},'Location','bestoutside')