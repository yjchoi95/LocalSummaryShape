function X = geodesic_sphere_Full(q1,q2,stp)


[N,T] = size(q1);
% Q = zeros(N,T,stp);
X = zeros(N,T,stp);
s = linspace(0,1,stp+1);%dist*(t-1)/6;

for t=1:(stp+1)
    X(:,:,t) = (1-s(t))*q1 + s(t)*q2;
    X(:,:,t) = ProjectC(X(:,:,t));
    %X(:,:,t) = q_to_curve(Q(:,:,t));
end
% theta = acos(InnerProd_Q(q1,q2));
% if theta > 0.0001
%     for t=1:stp+1
%         tau = (t-1)/stp;
%         X(:,:,t) = (sin((1-tau)*theta)*q1 + sin(tau*theta)*q2)/sin(theta);
%         X(:,:,t) = ProjectC_keepScale(X(:,:,t));      
%     end
% else
%     for t=1:stp+1
%         X(:,:,t) = q1;
%     end
% end

return;
