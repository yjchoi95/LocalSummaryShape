function beta_scaled = scaling(beta)

    dim = size(beta);
    T = dim(2);
    t = linspace(0,1,T);
    spacing = t(2)-t(1);
    beta_dot = gradient(beta, spacing);
    % compute the original length
    L = trapz(t,sqrt(sum(beta_dot.^2))); 
%     disp(strcat('original length: ', num2str(L)));
    beta_scaled = beta./L; %/10

%     % check if the length is 1
%     len = trapz(t,sqrt(sum(gradient(beta_scaled, spacing).^2)));
%     disp(strcat('after scaling: ', num2str(len)));