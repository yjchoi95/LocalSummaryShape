function fig = plot_local_summary_direction(ref, spatial_domain, xs, ys, ...
    psamp_mean, psamp_std, filename_to_save)


x_min = spatial_domain(1);
x_max = spatial_domain(2);
y_min = spatial_domain(3);
y_max = spatial_domain(4);

scale_shape = (x_max-x_min)/(length(unique(ref(:,1)))/1.8);
fig         =  figure; hold on;
scatter(ref(:,1),ref(:,2),'filled','MarkerEdgeColor','k','MarkerFaceColor','k');
for i = 1:size(psamp_mean, 3)
    s0x = ref(i, 1);
    s0y = ref(i, 2);
    
    % Scale psamp_mean for plotting
    scaled_mean_x = scale_shape * psamp_mean(1,:,i) + s0x;
    scaled_mean_y = scale_shape * psamp_mean(2,:,i) + s0y;
    
    % Plot mean curve
    plot(scaled_mean_x, scaled_mean_y, 'LineWidth', 1, 'Color', 'k');
    
    % Scale psamp_mean_std1 for plotting
    scaled_psamp_std1_x = scale_shape * psamp_std(1,:,i) + s0x;
    scaled_psamp_std1_y = scale_shape * psamp_std(2,:,i) + s0y;
    
    % obtain deviation
    scaled_dev_x = scaled_psamp_std1_x - scaled_mean_x;
    scaled_dev_y = scaled_psamp_std1_y - scaled_mean_y;
    
    % Plot deviation vectors as quivers
    idx = 1:3:length(scaled_mean_x);  % sample every 5rd index
    quiver(scaled_mean_x(idx), scaled_mean_y(idx), ...
       scaled_dev_x(idx), scaled_dev_y(idx), ...
       'Color', 'b', 'LineWidth', 1.5);

end

for i_x = 1:length(xs)
    line([xs(i_x) xs(i_x)], [y_min y_max], 'LineStyle', '--', 'Color', 'k');
end

for i_y = 1:length(ys)
    line([x_min x_max], [ys(i_y) ys(i_y)], 'LineStyle', '--', 'Color', 'k');
end
axis equal;
xlim([x_min,x_max]);
ylim([y_min,y_max]);

set(gca, 'FontSize', 25);
% xtick = get(gca, 'XTick');
% ytick = get(gca, 'YTick');
% set(gca, 'FontSize', 25);
% set(gca, 'XTick', xtick);
% set(gca, 'YTick', ytick);
hold off;
tightfig;

saveas(fig, strcat(filename_to_save,'.fig'));
saveas(fig, strcat(filename_to_save,'.png'));

