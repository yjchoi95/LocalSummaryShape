function [psamp,spat_info,spat_dist_vec] = ...
    GenerateCorrelatedSamples_Fourier(design,plotit,range)
% generate spatial shapes based on fourier model with mean zero and homogenous PP    
% Inputs: mu,X = outputs from FindMean3Dcurve, numsamp is number of
    
    
    % random samples to generate
    addpath('Functions');
    % folder0= 'Results';
    
    
    
    %%%%%%%%%%%%%%%%%%%% Parameter Setting %%%%%%%%%%%%%%%%%%%%
    sig_a = 0.7; % overall sd of matern covariance
    nu = 0.5; % smoothing parameter of matern covariance
    if nargin < 3
        range = 0.5; % range parameter of matern covariance
    end

    
    
    %%%%%%%%%%%%%%%%%%%% Spatial Design and PC coefficients %%%%%%%%%%%%%%%%%%%%
    if(design == 1)
        [x, y] = meshgrid(1:4); 
    else % design ==2
        x = 4*rand(50,1);
        y = 4*rand(50,1);
    end
    xy = [x(:), y(:)];
    spat_dist_vec = pdist(xy);
    N = size(xy, 1); % total number of observations
    spat_info = table((1:N)', xy(:,1), xy(:,2),'VariableNames', {'site', 'x', 'y'});

    
    rho = 4*range;%range*max(spat_dist_vec); 
    C_vec = MaternCovariance(spat_dist_vec, sig_a, nu, rho);
    C = squareform(C_vec) + diag(ones(N,1)*sig_a^2);

    
    
    %%%%%%%%%%%%%%%%%%%% Shape %%%%%%%%%%%%%%%%%%%%
    J = 4;
    T = 100;
    coef = mvnrnd(zeros(J*2,N),C);% + (mean_cluster).';
    %c = 4; % radius of the starting point
    psamp = zeros(2,T,N);
    
    theta = linspace(-pi, pi, T);
    f = @(x, coef, iter) coef(1)*cos(iter*x) + coef(2)*sin(iter*x);
    
    for i = 1:N
        r = zeros(T, J);
        coef1 = coef(:,i);
        for j = 1:J
            coef2 = coef1(2*(j-1)+1:2*(j-1)+2);
            r(:,j) = f(theta, coef2, j-1);
        end

        rsum = sum(r, 2);
        rsum = rsum + abs(min(rsum)) + abs(max(rsum));
        x = zeros(size(r,1),1);
        y = zeros(size(r,1),1);
        for j = 1:T
            x(j) = cos(theta(j))*rsum(j);
            y(j) = sin(theta(j))*rsum(j);
        end
        %     subplot(1, 2, 1);
        %     plot(theta, rsum); ylim([min(rsum)-10, max(rsum)+10]); grid on;
        %     xlabel('theta'); ylabel('rsum'); title('rsum vs theta');
        %     line([min(theta), max(theta)], [0, 0], 'Color', 'red');
        % 
        %     subplot(1, 2, 2);
        %     plot(x, y); hold on; grid on;
        %     plot(0, 0, 'ro', 'MarkerSize', 8); line([0, 0], ylim, 'Color', 'red');
        %     line(xlim, [0, 0], 'Color', 'red'); xlabel('x'); ylabel('y'); title('Circular trajectory');
        %     plot(x(1), y(1), 'bo', 'MarkerSize', 8);

        psamp(:,:,i) = [x,y].';
    end


    
    %%%%%%%%%%%%%%%%% Plotting spatial location %%%%%%%%%%%%%%%%%
    if(plotit == 1 && design == 1)
        colors = {'r', 'b', 'g', 'm'};
        plot_order = [4,8,12,16,3,7,11,15,2,6,10,14,1,5,9,13];

        figure;
        tiledlayout(4,4)
        for i = 1:N
            nexttile
            ind_col = spat_info.true_cluster(plot_order(i));
            plot(psamp(1,:,plot_order(i)),psamp(2,:,plot_order(i)), 'Color', colors{ind_col}, 'LineWidth',2); axis equal;
        end
    end

    if(plotit == 1 && design == 2)
        figure;
        scatter(spat_info.x,spat_info.y);
        % hold on;
        % xline(2, '--k', 'LineWidth', 2);
        % yline(2, '--k', 'LineWidth', 2);
        % legend('off')
        % legend('1', '2', '3', '4')
        % hold off; 
    end
end


