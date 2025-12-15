function [v,d,q2n] = ElasticShootingVector_rotation(q1,q2)

% reparamFlag=1 if using re-parameterization

q2n = Find_Seed_unique(q1,q2);
q2n = q2n/sqrt(InnerProd_Q(q2n,q2n)); % scaling

if InnerProd_Q(q1,q2n)>1 % handle imaginery part of d when the two shapes are equivalent
    d = 0;
else
    d = acos(InnerProd_Q(q1,q2n));
end

if d < 0.0001
    v = zeros(size(q1));
else
    v = (d/sin(d))*(q2n - cos(d)*q1);
end