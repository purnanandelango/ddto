function [] = plot_condensed_l1(sol,lin_mrk)

legend('AutoUpdate','on')
plot(sol.X{1}(1,:),sol.X{1}(2,:),lin_mrk{1},'DisplayName',lin_mrk{2},'Color',lin_mrk{3});
legend('AutoUpdate','off')
hold on
plot(sol.X{2}(1,:),sol.X{2}(2,:),lin_mrk{1},'Color',lin_mrk{3});
plot(sol.X{3}(1,:),sol.X{3}(2,:),lin_mrk{1},'Color',lin_mrk{3});
plot(sol.X{4}(1,:),sol.X{4}(2,:),lin_mrk{1},'Color',lin_mrk{3});

end