function [fig,min_dev,max_dev] = plot_local_summary_magnitude(ref, spatial_domain, xs, ys, ...
    psamp_mean, psamp_std, path_filename_to_save)


x_min = spatial_domain(1);
x_max = spatial_domain(2);
y_min = spatial_domain(3);
y_max = spatial_domain(4);


% obtain deviation and scale of deviation
psamp_dev = zeros(2,size(psamp_mean,2),size(ref,1));
for i = 1:size(ref,1)
    psamp_dev(:,:,i) = psamp_std(:,:,i) - psamp_mean(:,:,i);
end

scale_shape = (x_max-x_min)/(length(unique(ref(:,1)))/1.8);
fig         =  figure; hold on;
scatter(ref(:,1),ref(:,2),'filled','MarkerEdgeColor','k','MarkerFaceColor','k');
magnitude   = squeeze(sqrt(sum(psamp_dev.^2, 1))); % Euclidean norm of the vectors
max_dev     = max(magnitude,[],'all');
min_dev     = min(magnitude,[],'all');

for i = 1:size(psamp_mean, 3)
    s0x = ref(i,1);
    s0y = ref(i,2);
    fig = scatter(scale_shape * psamp_mean(1,:,i) + s0x, ...
            scale_shape * psamp_mean(2,:,i) + s0y, ...
            25, magnitude(:,i), 'filled');
end
colorbar;
colormap(jet);
clim([min_dev, max_dev]); 

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

saveas(fig, strcat(path_filename_to_save,'.fig'));
saveas(fig, strcat(path_filename_to_save,'.png'));


end