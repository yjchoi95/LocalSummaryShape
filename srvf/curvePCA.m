function [psamp,K] = curvePCA(q, muq, pc, path_filename)

b  = size(muq,2);
n  = size(q,3);
VV = zeros(n,b*2);
for i=1:n
    tmp = ElasticShootingVector(muq,q(:,:,i),1);
    VV(i,1:b) = tmp(1,:);
    VV(i,b+1:2*b) = tmp(2,:);
end

K = cov(VV);
[U,S,~] = svd(K);

T = length(VV(1,:))/2;
s = sqrt(diag(S));

k=1;
fig=figure; clf; hold on;
cols = [
    1 0 0;   % Red
    .9 0 .1;   % orange 
    .8 0 .2;   % green
    .7 0 .3;  % blue
    .6 0 .4;   % purple
    .5 0 .5;
    .4 0 .6
];
for t=-3:1:3
    u = t*s(pc)*U(:,pc)';
    UU(1,1:T) = u(1:T);
    UU(2,1:T) = u(T+1:2*T);
    UUn=UU;


    q = ElasticShooting(muq,UUn);
    p = ReSampleCurve(q_to_curve(q),T);
    psamp(:,:,k) = p - mean(p,2);
    
    plot(t * 0.45 + p(1,:), p(2,:),'color',cols(k,:),'LineWidth',3);
    axis equal off; view([1 90])
    saveas(fig, path_filename);
    k = k + 1;
end


end
% del = 0.3;
% mu = sum(q,3)/n;
% mu = mu/sqrt(InnerProd_Q(mu,mu));
% mu = ProjectC(mu);
% 
% for iter =1:Niter
%     disp(iter);
%     vm = 0;    
%     for i=1:n
%        % [iter i]
%         v = ElasticShootingVector(mu,q(:,:,i),1); % inv. of exp. mapping on mu_q for qi
%         vm = vm + v;
% 
%     end
%        vm = vm/n;
% 
%       % E(iter) = sqrt(InnerProd_Q(vm,vm));
%        E(iter) = norm(vm,'fro');
%        mu = ElasticShooting(mu,del*vm); % update mu 
% 
%        E
% 
% end
% 
% figure(101); 
% plot(E)
% 
% for i=1:n
%     [qbest,Rbest(:,:,i)] = Find_Rotation_and_Seed_unique(mu,q(:,:,i),1); % q with optimal rotation and reparam. to the mu updated above
%     qn(:,:,i) = qbest;
%     Xn(:,:,i) = q_to_curve(qbest);
% end
% % figure(21); clf; 
% muq = mu;
% mu = q_to_curve(mu);
% 
% p = mu;
% plot(p(1,:),p(2,:),'LineWidth',3);
% axis equal;

% reparamFlag=1;
% b=size(mu,2);
% 
% for i=1:n
%     i
%     tmp = ElasticShootingVector(mu,q(:,:,i),reparamFlag);
%     VV(i,1:b) = tmp(1,:);
%     VV(i,b+1:2*b) = tmp(2,:);
% end
% 
% K = cov(VV);
% [U,S,V] = svd(K);
% 
% T = length(VV(1,:))/2;
% s = sqrt(diag(S));
% 
% i=1;
% figure(11); clf; hold on;
% for t=-3:1:3
%     u = t*s(i)*U(:,i)';
%     UU(1,1:T) = u(1:T);
%     UU(2,1:T) = u(T+1:2*T);
%     UUn=UU;
% 
%     q = ElasticShooting(muq,UUn);
% 
% 
%     p = q_to_curve(q);
% 
% 
%     if (t==0)
%         plot(t*0.32+p(1,:),p(2,:),'r','LineWidth',3);
%     else
%         plot(t*0.32+p(1,:),p(2,:),'LineWidth',3);
%     end
%     axis equal off; view([1 90])
% end
% 
% i = 2;
% figure(12); clf; hold on;
% for t=-3:1:3
%     u = t*s(i)*U(:,i)'/2;
%     UU(1,1:T) = u(1:T);
%     UU(2,1:T) = u(T+1:2*T);
%     UUn=UU;
% 
%     q = ElasticShooting(muq,UUn);   
%     p = q_to_curve(q);
% 
% 
%     if (t==0)
%         plot(t*0.31+p(1,:),p(2,:),'r','LineWidth',3);
%     else
%         plot(t*0.31+p(1,:),p(2,:),'LineWidth',3);
%     end
%     axis equal off; view([1 90])
% end

% i = 3;
% figure(13); clf; hold on;
% for t=-3:1:3
%     u = t*s(i)*U(:,i)'/2;
%     UU(1,1:T) = u(1:T);
%     UU(2,1:T) = u(T+1:2*T);
%     UUn=UU;
%     
%     q = ElasticShooting(mu,UUn);   
%     p = q_to_curve(q);
%     
%     if (t==0)
%         plot(t*0.31+p(1,:),p(2,:),'r','LineWidth',3);
%     else
%         plot(t*0.31+p(1,:),p(2,:),'LineWidth',3);
%     end
%     axis equal off; view([1 90])
% end
% 
% mun=mean(VV);
% idx=1;
% figure(17); clf; hold on;
% while(idx<7)
%         v = mun + zeros(1,2*T);
%         for i=1:b
%             % Gaussian
%             c(i) = randn*s(i);
%             v = v + c(i)*U(:,i)';
%         end
%         
%         vn(1,1:T) = v(1:T);
%         vn(2,1:T) = v(T+1:2*T);
%         qsamp = ElasticShooting(mu,vn);
%         psamp = q_to_curve(qsamp);
% 
%         z = plot(psamp(1,:)+.31*idx,psamp(2,:),'LineWidth',3);
%         axis equal off; view([1 90]);
%         idx=idx+1;
% end