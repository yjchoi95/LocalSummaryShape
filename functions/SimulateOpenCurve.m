function [data,coord,distance] = SimulateOpenCurve(design,plotit,range)

    folder0='Results';
    addpath('Functions');
    
    %%%%%%%%%%%%%%%%%%%%% parameter definition %%%%%%%%%%%%%%%%%%%%%
    N_t = 50; % number of points in each shape
    t = linspace(0, 1, N_t); % timestamp
    sig_a = 1; % overall sd of matern covariance
    sig_e = .5; % sd of gaussian noise
    if nargin < 3
        range = 0.5; % range parameter of matern covariance
    end
    scale_b = 1; % range parameter of correlated uniform distribution [-B,B]
    % peak and skewness parameter
    peak_diff = 2;
    skew_diff = .8;
    
    mu = betapdf(t,8,2);
    %mu = normpdf(t, 0.5, 0.1); % true mean shape
     %mu = -cos(t*(2*pi));
    % plot(t, mu, 'LineWidth', 1);

    if(design == 1)
        [x, y] = meshgrid(1:4); 
    else % design ==2
        x = 4*rand(30,1);
        y = 4*rand(30,1);
    end
    coordinates = [x(:), y(:)];
    distance = pdist(coordinates);
    N = size(coordinates, 1); % total number of observations
    coord = table((1:N)', coordinates(:,1), coordinates(:,2),'VariableNames', {'site', 'x', 'y'});

    % design. Define clusters
    peak_cluster = zeros(N,1);
    peak_cluster(coord.x>2 & coord.y>2) = 4;
    peak_cluster(coord.x>2 & coord.y<=2) = 3;
    peak_cluster(coord.x<=2 & coord.y>2) = 2;
    peak_cluster(coord.x<=2 & coord.y<=2) = 1;

    % coord.peak_cluster = categorical(peak_cluster);
    coord.peak_cluster = peak_cluster;
    coord.skewness_cluster = coord.peak_cluster;
        
    if(plotit == 1)
        figure;
        gscatter(coord.x,coord.y,coord.peak_cluster);
        hold on;
        if(design == 1)
            bound=2.5;
        else
            bound=2;
        end
        xline(bound, '--k', 'LineWidth', 2);
        yline(bound, '--k', 'LineWidth', 2);
        legend('off')
        legend('1', '2', '3', '4')
        hold off; 
    end

    % compute matern covariance: sig_a-var, rho-range
    v = 0.5;
    rho = range*max(distance);
    C_vec = MaternCovariance(distance, sig_a, v, rho);
    C = squareform(C_vec) + diag(ones(N,1)*sig_a^2);

    %%%%%%%%%%%%%%%%%%%%% differ peak %%%%%%%%%%%%%%%%%%%%%
    a = mvnrnd(repmat(5,1,N),C) + (coord.peak_cluster).'*peak_diff;%double(coord.peak_cluster).'*peak_diff;
    curve = mu.'*a;

    %%%%%%%%%%%%%%%%%%%%% differ skewness %%%%%%%%%%%%%%%%%%%%%

    % calculate beta_b
    beta_b = exp(normcdf(mvnrnd(zeros(1,N),C./sig_a^2),0,1).*2*scale_b - scale_b + coord.skewness_cluster.'.*skew_diff);

    % get gamma functions
    gam = zeros(N_t,N);
    for i = 1:N
      gam(:,i) = betacdf((t-min(t))/(max(t)-min(t)),1,beta_b(i));
    end

%     t_gam = zeros(N_t,N);
    % reparameterize the curves
    for i = 1:N
    %     figure;
    %     plot(t,curve(:,i));
    %     ylim([0,60]);
    %     
    %     figure;
    %     plot(t,gam(:,i));
    %     
        curve(:,i) = Group_Action_by_Gamma_Coord(curve(:,i).',gam(:,i));

    %     figure;
    %     plot(t,curve(:,i));
    %     ylim([0,60]);
    end

    % add some gaussian noise
    normal_error = normrnd(0,sig_e,[N_t,N]);
    data_matrix = curve + normal_error;
    %data_matrix = curve;
    % plot the curves accounting for spatial info
    if(plotit == 1 && design == 1)
        colors = {'r', 'b', 'g', 'm'};
        plot_order = [4,8,12,16,3,7,11,15,2,6,10,14,1,5,9,13];

        figure;
        tiledlayout(4,4)
        for i = 1:N
            nexttile
            ind_col = coord.peak_cluster(plot_order(i));
            plot(t,data_matrix(:,plot_order(i)), 'Color', colors{ind_col});
            ylim([0,50]);
        end
    end

    data = zeros(2,N_t,N);
    for i = 1:N
        data(1,:,i) = t;
        data(2,:,i) = data_matrix(:,i);
    end

    save(strcat(folder0,'/simulated_data.mat'),'data');
    save(strcat(folder0,'/spatial_info.mat'),'coord');
    save(strcat(folder0,'/spatial_distance.mat'),'distance');

end