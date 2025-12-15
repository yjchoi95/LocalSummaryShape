function [q] = curve_to_q(p)

[n,N] = size(p);
for i = 1:n
    v(i,:) = gradient(p(i,:),1/(N)); % v = beta dot at t=i
end

for i = 1:N
    L(i) = sqrt(norm(v(:,i),'fro')); % L = sqrt of speed at t=i 
    if L(i) > 0.00001
        q(:,i) = v(:,i)/L(i); % srvf representation
    else
        q(:,i) = 0*ones(n,1);
    end
end

q = q/sqrt(InnerProd_Q(q,q)); % scale q to have a length of 1