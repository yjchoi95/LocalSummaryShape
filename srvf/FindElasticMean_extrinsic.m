function [muX,muq,Xn,qn,delta_norm,Rbest,muq_seq,q] = FindElasticMean_extrinsic(data_matrix)

stp     = 3;
eps     = 5e-4;
maxIter = 40;
[~,T,n] = size(data_matrix);
q = zeros(size(data_matrix));
for i=1:n
    X = ReSampleCurve(data_matrix(:,:,i),T);
    q(:,:,i) = curve_to_q(X);
end

muq_seq = mean(q,3);
muq_seq = muq_seq/sqrt(InnerProd_Q(muq_seq,muq_seq));
iter    = 1;
qn    = zeros(size(q));
Rbest = zeros(2,2,n);
while iter<=maxIter

    disp(iter);
    q_current = muq_seq(:,:,iter);
    % alignment
    for i=1:n
        [q_aligned,Rbest(:,:,i)] = ...
            Find_Rotation_and_Seed_unique(q_current,q(:,:,i),1);
        qn(:,:,i) = q_aligned/sqrt(InnerProd_Q(q_aligned,q_aligned)); % scaling
    end
    
    q_current_rep     = repmat(q_current, [1, 1, size(q_aligned, 3)]);  % q_current is 2x100
    deviations        = q_aligned - q_current_rep;
    delta             = sum(deviations, 3)/n;  % same shape as q_current
    z                 = q_current + stp * delta;
    q_new             = ProjectC(z); 
    muq_seq(:,:,iter+1) = q_new;
    delta_norm(iter)  = sqrt(InnerProd_Q(delta,delta));
    disp(delta_norm);
    if delta_norm(iter) < eps
        break
    end
    iter = iter+1;
end

Xn = zeros(size(q));
for i=1:n
    qbest     = Find_Rotation_and_Seed_unique(q_new,q(:,:,i),1);
    qn(:,:,i) = qbest/sqrt(InnerProd_Q(qbest,qbest)); % scaling
    Xn(:,:,i) = q_to_curve(qbest);
end
muq = q_new;
muX = q_to_curve(muq);

end