function gamI = invertGamma(gam)

N = length(gam);
x = [1:N]/N;
% xx = [1:N-1]/(N-1);
% gamI1 = spline(gam,x,x);
gamI = interp1(gam,x,x,'linear');

% gamI = gamI/gamI(N);
