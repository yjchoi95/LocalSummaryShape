function [muX,muq,Xn,qn,E,Rbest,muq_seq] = FindElasticMean_extrinsic0(data_matrix)


eps = 1e-4;
maxIter = 20;
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
    
    % extrinsic mean
    z = mean(qn,3);
    
    % put it back to the sphere space
    q_new = ProjectC(z); % it scales the q 
    muq_seq(:,:,iter+1) = q_new;
    E(iter) = sqrt(InnerProd_Q(q_new-q_current, q_new-q_current));
    disp(E);
    if E(iter) < eps
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