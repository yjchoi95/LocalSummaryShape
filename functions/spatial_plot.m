function fig=spatial_plot(data_matrix, xy, true_cluster)
    
    % scale the length to 1
    data_matrix_length1 = scaling(data_matrix).*1;
    % data_matrix_length1=data_matrix;
    
    % find centroid
    centroid = Find_Centroids(data_matrix_length1, 0);
    
    % center it to zero
    centroid_reshaped = reshape(centroid.', [2, 1, size(centroid, 1)]);
    center_at_zero = data_matrix_length1 - centroid_reshaped;
    
    % relocate to their spatial location
    coord_reshaped = reshape(xy.', [2, 1, size(xy, 1)]);
    center_adj = center_at_zero + coord_reshaped;
    
    % plot
    colors = remove_yellow_from_jet(length(unique(true_cluster)));
    
    fig=figure; hold on;
    for i = 1:size(data_matrix,3)
        label = true_cluster(i);
        plot(center_adj(1,:,i), center_adj(2,:,i), 'Linewidth', 1.5, 'Color', colors(label,:)); 
        plot(xy(i,1), xy(i,2), 'k.', 'Linewidth', 3);
    end
    axis equal;
    hold off;
    %saveas(gcf,'psamp.png');
    
    
    function modified_jet = remove_yellow_from_jet(n)
        % Get the jet colormap
        J = jet(n);
    
        % Identify the yellow region (yellow has high red and green values)
        yellow_region = (J(:,1) > 0.5 & J(:,2) > 0.5 & J(:,3) < 0.5);
    
        % Remove the yellow region
        J(yellow_region, :) = [];
    
        %   Interpolate to get the colormap of size n
        xi = linspace(1, length(J), n);
        modified_jet = interp1(1:length(J), J, xi);
    end
end


% % example
% [data_matrix, coord, ~] = GenerateCorrelatedSamples(2,0,'bone',1);
% coordinates = [coord.x, coord.y];
% true_cluster = coord.true_cluster;
% spatial_plot(data_matrix, coordinates, true_cluster);