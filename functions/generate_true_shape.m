function [psamp] = generate_true_shape(coef)

    T = 100; % sample points       
    J = length(coef)/2;
    r = zeros(T, J); 
    theta     = linspace(-pi, pi, T);
    f         = @(x, coef, iter) coef(1)*cos(iter*x) + coef(2)*sin(iter*x);

    for j = 1:J
        coef0  = coef(2*(j-1)+1:2*(j-1)+2);
        r(:,j) = f(theta, coef0, j-1);
    end

    rsum      = sum(r, 2);
    rsum      = rsum + abs(min(rsum)) + abs(max(rsum));

    x = zeros(size(r,1),1);
    y = zeros(size(r,1),1);
    for j = 1:T
        x(j) = cos(theta(j))*rsum(j);
        y(j) = sin(theta(j))*rsum(j);
    end
    psamp = [x,y].';             
    psamp = scaling(psamp);
end


