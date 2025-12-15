function [v,d,q2n] = ElasticShootingVector(q1,q2,reparamFlag)

% reparamFlag=1 if using re-parameterization

q2n = Find_Rotation_and_Seed_unique(q1,q2,reparamFlag);
q2n = q2n/sqrt(InnerProd_Q(q2n,q2n)); % scaling

d = acos(InnerProd_Q(q1,q2n));
if d < 0.0001
    v = zeros(size(q1));
else
    v = (d/sin(d))*(q2n - cos(d)*q1);
end