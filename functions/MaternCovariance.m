% D: pdist object
% sigma: sqrt of variance
% nu: smoothing parameter, usually =0.5
% range: range parameter
% function C = MaternCovariance(D, sigma, nu, range)
%     C = sigma^2 * 1/(2^(nu-1) * gamma(nu)) * (D/range).^nu.*besselk(nu,D/range);
% end


function val = MaternCovariance(D, sigma, nu, range)
    val = zeros(size(D));
    % on the diagonal
    val(D == 0) = sigma^2;
    
    % covariance for nonzero distances
    idx = D > 0;
    factor = sqrt(2 * nu) * D(idx) / range;
    val(idx) = sigma^2 * (1 / (2^(nu - 1) * gamma(nu))) * ...
                   (factor.^nu) .* besselk(nu, factor);

end