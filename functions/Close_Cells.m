function [closed_cells] = Close_Cells(unclosed_cells, plotit)
% input: unclosed_cells - d by T by N cell data after ReSampleCurve
% input: plotit - plot or not
% output: closed_cells - d by T by N cell data whose curves are closed
%plotit = true;

[d,T,N] = size(unclosed_cells);
closed_cells = zeros(d,T,N);

for i = 1:N
    % close curve
    q = curve_to_q_unscaled(unclosed_cells(:,:,i));
    q_new = ProjectC_unscaled(q);
    p_new = q_to_curve(q_new);
    
    % obtain spatial difference between the original curve and closed curve
    x_diff = min(unclosed_cells(1,:,i)) - min(p_new(1,:));
    y_diff = min(unclosed_cells(2,:,i)) - min(p_new(2,:));
    
    % relocate closed curve to the original spatial location
    % w.r.t. min of x and y
    closed_cells(1,:,i) = p_new(1,:) + x_diff;
    closed_cells(2,:,i) = p_new(2,:) + y_diff;
 
end

% plot if plotit=true;
if(plotit)
    figure; 
    plot(closed_cells(1,:,1), closed_cells(2,:,1));
    xlim([-100,1100]); ylim([-100,1100]);
    hold on;
    for i = 2:N
        plot(closed_cells(1,:,i), closed_cells(2,:,i));
    end
    title('Cells - Closed Curves');
    hold off;
    
    figure;
    plot(unclosed_cells(1,:,1), unclosed_cells(2,:,1));
    hold on;
    for i = 2:N
        plot(unclosed_cells(1,:,i), unclosed_cells(2,:,i));
    end
    title('Cells - Original');
    xlim([-100,1100]); ylim([-100,1100]);
    hold off;
end

end
