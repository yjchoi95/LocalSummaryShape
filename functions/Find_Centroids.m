function s = Find_Centroids(closed_cells, plotit)
% input: d by T by N closed cells, where d=dimension, T= sample points, 
% and N = number of curves
% output: centroids of curves 

    % allocate space for centroids
    s = zeros(size(closed_cells,3),2);

    for i = 1:size(closed_cells,3)
    % find area enclosed by a simple curve 
    curve = closed_cells(:,:,i);

    x = curve(1,:);
    y = curve(2,:);

    % ensure the vectors are column vectors
    if isrow(x), x = x'; end
    if isrow(y), y = y'; end

    % append the first point to the end to close the polygon
    x = [x; x(1)];
    y = [y; y(1)];

    % number of points
    n = length(x);

    % compute the area
    A = 0.5 * (sum(x(1:end-1).*y(2:end) - x(2:end).*y(1:end-1)));

    % compute the centroids
    x_bar = (1/(6*A)) * sum((x(1:end-1) + x(2:end)) .* (x(1:end-1).*y(2:end) - x(2:end).*y(1:end-1)));
    y_bar = (1/(6*A)) * sum((y(1:end-1) + y(2:end)) .* (x(1:end-1).*y(2:end) - x(2:end).*y(1:end-1)));

    s(i,:) = [(x_bar),(y_bar)];
    end

    % plotting
    if(plotit)
        figure; 
        plot(closed_cells(1,:,1), closed_cells(2,:,1));
        xlim([-100,1100]); ylim([-100,1100]);
        hold on;
        plot(s(1,1), s(1,2), '.');
        for i = 2:size(closed_cells,3)
            plot(closed_cells(1,:,i), closed_cells(2,:,i));
            plot(s(i,1), s(i,2), '.');
        end
        title('Cells - Closed Curves with centroids');
        hold off;
    end

end