function centers = compute_centers(subregions)
    % compute_centers calculates the center of each subregion.
    % Input:
    %   subregions - a matrix where each row is [x_min, x_max, y_min, y_max]
    % Output:
    %   centers - a matrix where each row is [center_x, center_y] for the subregion

    % Vectorized computation:
    center_x = (subregions(:,1) + subregions(:,2)) / 2;
    center_y = (subregions(:,3) + subregions(:,4)) / 2;
    
    centers = [center_x, center_y];
end
