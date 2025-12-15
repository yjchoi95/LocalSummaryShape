function [fig,center_adj]=spatial_plot_shape(data_matrix, xy)
        
    % find centroid
    centroid = Find_Centroids(data_matrix, 0);
    
    % center it to zero
    centroid_reshaped = reshape(centroid.', [2, 1, size(centroid, 1)]);
    center_at_zero = data_matrix - centroid_reshaped;
    
    % relocate to their spatial locations
    coord_reshaped = reshape(xy.', [2, 1, size(xy, 1)]);
    center_adj = center_at_zero + coord_reshaped;
    
    x_min = min(min(center_adj(1,:,:))); 
    x_max = max(max(center_adj(1,:,:))); 
    y_min = min(min(center_adj(2,:,:))); 
    y_max = max(max(center_adj(2,:,:)));

    fig=figure; hold on;
    for i = 1:size(data_matrix,3)
        plot(center_adj(1,:,i), center_adj(2,:,i), 'k', 'Linewidth', 1.5); 
        plot(xy(i,1), xy(i,2), 'k.', 'Linewidth', 3);
    end
    axis equal;
    xlim([x_min,x_max]);
    ylim([y_min,y_max]);
    set(gca,'FontSize',20);
    % tightfig;
    hold off;
    
end

