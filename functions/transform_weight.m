function w = transform_weight(lambda)
    abs_min = abs(min(lambda));
    denom = sum(lambda+abs_min);
    w = (lambda+abs_min)./denom;
end