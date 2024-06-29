function [] = plot_condensed(p1,sol,p2,sol2,p3,sol3,trgt_ord)

mrk_scl = 20;
anim_rate = 0.2; % 1/fps
% arrow_scl = 0.15;
% stem_width = 3;
file_name = 'traj_anim_3.gif';
save_anim = true;
img = [];

interpreter = "tex";
% interpreter = "latex";

fig = figure('Color',[1 1 1]);

% initpt = plot3(p1.z0(1),p1.z0(2),p1.z0(3),'+k','MarkerSize',mrk_scl);
hold on
view(91,17);
% trgt1 = plot3(p1.zf(1,1),p1.zf(2,1),p1.zf(3,1),'dk','MarkerSize',mrk_scl);
% trgt2 = plot3(p1.zf(1,2),p1.zf(2,2),p1.zf(3,2),'ok','MarkerSize',mrk_scl);
% trgt3 = plot3(p1.zf(1,3),p1.zf(2,3),p1.zf(3,3),'<k','MarkerSize',mrk_scl);
% trgt4 = plot3(p1.zf(1,4),p1.zf(2,4),p1.zf(3,4),'sk','MarkerSize',mrk_scl);

% legend([initpt,trgt1,trgt2,trgt3,trgt4],{'Source','Target 1','Target 2','Target 3','Target 4'},...
%        'AutoUpdate','off',...
%        'Position',[0.75, 0.50, .1, .1]);
ax = gca;
grid on
ax.DataAspectRatioMode = 'manual';
ax.DataAspectRatio = [1 1 1];
ax.PlotBoxAspectRatioMode = 'manual';
ax.PlotBoxAspectRatio = [1 1 1];
ax.XLim = [1.1*min(min(p1.zf(1:2,:))) 1.1*max(max(p1.zf(1:2,:)))];
ax.YLim = [1.1*min(min(p1.zf(1:2,:))) max(max(p1.zf(1:2,:)))];
ax.ZLim = [-1 max(max(p1.z0(3,:)))];

xlabel('[m]');
ylabel('[m]');
zlabel('[m]');

% 1st deferrable segment
for j = 1:sol.dd_idx+1
    plot3(sol.X{1}(1,1:j),...
          sol.X{1}(2,1:j),...
          sol.X{1}(3,1:j),...
          'Color','black','Marker','.','MarkerSize',mrk_scl);

    % if j<=sol.dd_idx
    %     q = quiver3(sol.X{1}(1,j),...
    %                 sol.X{1}(2,j),...
    %                 sol.X{1}(3,j),...
    %                 sol.U{1}(1,j),...
    %                 sol.U{1}(2,j),...
    %                 sol.U{1}(3,j),arrow_scl);
    %     q.ShowArrowHead = 'off';
    %     q.Marker = '.';
    %     q.LineWidth = stem_width;
    %     q.Color = [0.6,0.5,0.5];
    % end 

    drawnow;
    frame = getframe(fig);
    img{end+1} = frame2im(frame);
end
idx = find(p1.tag==trgt_ord(1));
% 1st branch point
% plot3(sol.X{idx}(1,sol.dd_idx+1),...
%           sol.X{idx}(2,sol.dd_idx+1),...
%           sol.X{idx}(3,sol.dd_idx+1),...
%           'LineWidth',1.5,'Color','black','Marker','x');
[~,N] = size(sol.X{idx});
% Trajectory from 1st branch point to 1st rejected target
for j = sol.dd_idx+1:N
    plot3(sol.X{idx}(1,sol.dd_idx+1:j),...
          sol.X{idx}(2,sol.dd_idx+1:j),...
          sol.X{idx}(3,sol.dd_idx+1:j),...
          'Color','blue','Marker','.','MarkerSize',mrk_scl);

    % if j<=N-1
    %     q = quiver3(sol.X{idx}(1,j),...
    %                 sol.X{idx}(2,j),...
    %                 sol.X{idx}(3,j),...
    %                 sol.U{idx}(1,j),...
    %                 sol.U{idx}(2,j),...
    %                 sol.U{idx}(3,j),arrow_scl);
    %     q.ShowArrowHead = 'off';
    %     q.Marker = '.';
    %     q.LineWidth = stem_width;
    %     q.Color = [0 0.25 0.75];
    % end 

    drawnow;
    frame = getframe(fig);
    img{end+1} = frame2im(frame);    
end

% 2nd deferrable segment
for j = 1:sol2.dd_idx+1
    plot3(sol2.X{1}(1,1:j),...
          sol2.X{1}(2,1:j),...
          sol2.X{1}(3,1:j),...
          'LineWidth',1.5,'Color','black','Marker','.','MarkerSize',mrk_scl);
    % if j<=sol2.dd_idx
    %     q = quiver3(sol2.X{1}(1,j),...
    %                 sol2.X{1}(2,j),...
    %                 sol2.X{1}(3,j),...
    %                 sol2.U{1}(1,j),...
    %                 sol2.U{1}(2,j),...
    %                 sol2.U{1}(3,j),arrow_scl);
    %     q.ShowArrowHead = 'off';
    %     q.Marker = '.';
    %     q.LineWidth = stem_width;
    %     q.Color = [0.6,0.5,0.5];
    % end

    drawnow;
    frame = getframe(fig);
    img{end+1} = frame2im(frame);
end
idx = find(p2.tag==trgt_ord(2));
% 2nd branch point
plot3(sol2.X{idx}(1,sol2.dd_idx+1),...
      sol2.X{idx}(2,sol2.dd_idx+1),...
      sol2.X{idx}(3,sol2.dd_idx+1),...
      'LineWidth',1.5,'Color','black','Marker','x');
[~,N] = size(sol2.X{idx});
% Trajectory from 2nd branch point to 2nd rejected target
for j = sol2.dd_idx+1:N
    plot3(sol2.X{idx}(1,sol2.dd_idx+1:j),...
          sol2.X{idx}(2,sol2.dd_idx+1:j),...
          sol2.X{idx}(3,sol2.dd_idx+1:j),...
          'LineWidth',2,'Color','red','Marker','.','MarkerSize',mrk_scl);

    % if j<=N-1
    %     q = quiver3(sol2.X{idx}(1,j),...
    %                 sol2.X{idx}(2,j),...
    %                 sol2.X{idx}(3,j),...
    %                 sol2.U{idx}(1,j),...
    %                 sol2.U{idx}(2,j),...
    %                 sol2.U{idx}(3,j),arrow_scl);
    %     q.ShowArrowHead = 'off';
    %     q.Marker = '.';
    %     q.LineWidth = stem_width;
    %     q.Color = [0.75 0 0.25];
    % end

    drawnow;
    frame = getframe(fig);
    img{end+1} = frame2im(frame);
end

% 3rd deferrable segment
for j = 1:sol3.dd_idx+1
    plot3(sol3.X{1}(1,1:j),...
          sol3.X{1}(2,1:j),...
          sol3.X{1}(3,1:j),...
          'LineWidth',1.5,'Color','black','Marker','.','MarkerSize',mrk_scl);

    % if j<=sol3.dd_idx
    %     q = quiver3(sol3.X{1}(1,j),...
    %                 sol3.X{1}(2,j),...
    %                 sol3.X{1}(3,j),...
    %                 sol3.U{1}(1,j),...
    %                 sol3.U{1}(2,j),...
    %                 sol3.U{1}(3,j),arrow_scl);
    %     q.ShowArrowHead = 'off';
    %     q.Marker = '.';
    %     q.LineWidth = stem_width;
    %     q.Color = [0.6,0.5,0.5];
    % end    

    drawnow;
    frame = getframe(fig);
    img{end+1} = frame2im(frame);
end
idx = find(p3.tag==trgt_ord(3));
% 3rd branch point
plot3(sol3.X{idx}(1,sol3.dd_idx+1),...
      sol3.X{idx}(2,sol3.dd_idx+1),...
      sol3.X{idx}(3,sol3.dd_idx+1),...
      'LineWidth',1.5,'Color','black','Marker','x');
[~,N] = size(sol3.X{idx});
% Trajectory from 3rd branch point to 3rd rejected target
for j = sol3.dd_idx+1:N
    plot3(sol3.X{idx}(1,sol3.dd_idx+1:j),...
          sol3.X{idx}(2,sol3.dd_idx+1:j),...
          sol3.X{idx}(3,sol3.dd_idx+1:j),...
          'LineWidth',2,'Color','green','Marker','.','MarkerSize',mrk_scl);

    % if j<=N-1
    %     q = quiver3(sol3.X{idx}(1,j),...
    %                 sol3.X{idx}(2,j),...
    %                 sol3.X{idx}(3,j),...
    %                 sol3.U{idx}(1,j),...
    %                 sol3.U{idx}(2,j),...
    %                 sol3.U{idx}(3,j),arrow_scl);
    %     q.ShowArrowHead = 'off';
    %     q.Marker = '.';
    %     q.LineWidth = stem_width;
    %     q.Color = [0.25 0.75 0];
    % end    

    drawnow;
    frame = getframe(fig);
    img{end+1} = frame2im(frame);
end 

idx = find(p3.tag==trgt_ord(4));
[~,N] = size(sol3.X{idx});
% Trajectory from 3rd branch point to 4th rejected target
for j = sol3.dd_idx+1:N
    plot3(sol3.X{idx}(1,sol3.dd_idx+1:j),...
          sol3.X{idx}(2,sol3.dd_idx+1:j),...
          sol3.X{idx}(3,sol3.dd_idx+1:j),...
          'LineWidth',2,'Color','magenta','Marker','.','MarkerSize',mrk_scl);

    % if j<=N-1
    %     q = quiver3(sol3.X{idx}(1,j),...
    %                 sol3.X{idx}(2,j),...
    %                 sol3.X{idx}(3,j),...
    %                 sol3.U{idx}(1,j),...
    %                 sol3.U{idx}(2,j),...
    %                 sol3.U{idx}(3,j),arrow_scl);
    %     q.ShowArrowHead = 'off';
    %     q.Marker = '.';
    %     q.LineWidth = stem_width;
    %     q.Color = [0.75 0.25 0.75];
    % end

    drawnow;
    frame = getframe(fig);
    img{end+1} = frame2im(frame);
end

if save_anim
    for idx = 1:length(img)
        [A,map] = rgb2ind(img{idx},256);
        if idx == 1
            imwrite(A,map,file_name,'gif','LoopCount',Inf,'DelayTime',anim_rate);
        else
            imwrite(A,map,file_name,'gif','WriteMode','append','DelayTime',anim_rate);
        end
    end
end
% gif2avi(file_name,'.mp4')

end

%%% COMET CODE
% comet3(sol.X{1}(1,1:sol.dd_idx+1),sol.X{1}(2,1:sol.dd_idx+1),sol.X{1}(3,1:sol.dd_idx+1));
% plot3(p1.z0(1),p1.z0(2),p1.z0(3),'LineWidth',1.5,'Color','black','Marker','x');
% idx = find(p1.tag==trgt_ord(1));
% comet3(sol.X{idx}(1,sol.dd_idx+1:end),sol.X{idx}(2,sol.dd_idx+1:end),sol.X{idx}(3,sol.dd_idx+1:end));
% 
% comet3(sol2.X{1}(1,1:sol2.dd_idx+1),sol2.X{1}(2,1:sol2.dd_idx+1),sol2.X{1}(3,1:sol2.dd_idx+1));
% plot3(p2.z0(1),p2.z0(2),p2.z0(3),'LineWidth',1.5,'Color','black','Marker','x');
% idx = find(p2.tag==trgt_ord(2));
% comet3(sol2.X{idx}(1,sol2.dd_idx+1:end),sol2.X{idx}(2,sol2.dd_idx+1:end),sol2.X{idx}(3,sol2.dd_idx+1:end));
% 
% comet3(sol3.X{1}(1,1:sol3.dd_idx+1),sol3.X{1}(2,1:sol3.dd_idx+1),sol3.X{1}(3,1:sol3.dd_idx+1));
% plot3(p3.z0(1),p3.z0(2),p3.z0(3),'LineWidth',1.5,'Color','black','Marker','x');
% idx = find(p3.tag==trgt_ord(3));
% comet3(sol3.X{idx}(1,sol3.dd_idx+1:end),sol3.X{idx}(2,sol3.dd_idx+1:end),sol3.X{idx}(3,sol3.dd_idx+1:end));
% 
% comet3(sol3.X{idx}(1,sol3.dd_idx+1),sol3.X{idx}(2,sol3.dd_idx+1),sol3.X{idx}(3,sol3.dd_idx+1));
% idx = find(p3.tag==trgt_ord(4));
% comet3(sol3.X{idx}(1,sol3.dd_idx+1:end),sol3.X{idx}(2,sol3.dd_idx+1:end),sol3.X{idx}(3,sol3.dd_idx+1:end));
