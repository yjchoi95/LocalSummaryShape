function [psamp,thetas,scaled_biv_normal] = spatial_shape(property,dependency,xy,sig,l,range)
    % property: a struct with rotation, and scale
    % dependency: a struct with shape, rotation, and scale
    % xy: N by 2 spatial locations
    % range: range parameter of matern covariance
    % psamp: d by N
    
    N        = size(xy, 1); % total number of observations
    % matern covariance and generate coefficients
    D = pdist2(xy,xy); 
    C = MaternCovariance(D, sig, l, range);
    
    J = 3;
    T = 100; % sample points   
    if dependency.shape == 1
        coef = mvnrnd(zeros(J*2,N),C);
    else
        cov_mat          = eye(N) * sig^2;
        cov_coef         = .7;
        cov_mat(~eye(N)) = sig^2 * cov_coef;
        coef             = mvnrnd(zeros(J*2,N),cov_mat);
    end 
    psamp     = zeros(2,T,N);
    rsum      = zeros(T,N);
    rsum_orig = zeros(T,N);
    theta     = linspace(-pi, pi, T);
    f         = @(x, coef, iter) coef(1)*cos(iter*x) + coef(2)*sin(iter*x);

    for i = 1:N
        r = zeros(T, J);
        coef1 = coef(:,i);
        for j = 1:J
            coef2  = coef1(2*(j-1)+1:2*(j-1)+2);
            r(:,j) = f(theta, coef2, j-1);
        end

        rsum(:,i)      = sum(r, 2);
        rsum_orig(:,i) = rsum(:,i);
        rsum(:,i)      = rsum(:,i) + abs(min(rsum(:,i))) + abs(max(rsum(:,i)));

        x = zeros(size(r,1),1);
        y = zeros(size(r,1),1);
        for j = 1:T
            x(j) = cos(theta(j))*rsum(j,i);
            y(j) = sin(theta(j))*rsum(j,i);
        end

        psamp(:,:,i) = [x,y].';             
    end
    
    for j = 1:size(psamp,3)
        psamp(:,:,j) = scaling(psamp(:,:,j));
    end
   
    scaled_biv_normal = zeros(N,2);
    % rotation
    if property.rotation == 1
        if dependency.rotation == 1       
            biv_normal = mvnrnd(zeros(1,N),C,2);
            thetas     = zeros(1,N);
            
            for j = 1:N
               v                      = biv_normal(:,j)./norm(biv_normal(:,j),'fro');
               scaled_biv_normal(j,:) = v;
               angle                  = atan2(v(2), v(1));
               % Adjust angle to be in the range [0, 2*pi] if negative
               if angle < 0
                    angle = angle + 2*pi;
               end
               thetas(j) = angle;
            end

        else
            thetas = normcdf(mvnrnd(zeros(1, N), eye(N)))*2*pi;
        end

        for i = 1:size(psamp,3)
            theta        = thetas(i);
            O            = [cos(theta), -sin(theta); sin(theta), cos(theta)];
            psamp(:,:,i) = O*psamp(:,:,i);
        end
    end
    % scaling
    if property.scale == 1
        if dependency.scale == 1
            scales = normcdf(mvnrnd(zeros(1, N), C)) * .4 + .2;
        else
            scales = normcdf(mvnrnd(zeros(1, N), eye(N))) * .4 + .2;
        end

        for i = 1:size(psamp,3)
            scaled_psamp = scaling(psamp(:,:,i));
            psamp(:,:,i) = scales(i).*scaled_psamp;
        end        
    end
end


