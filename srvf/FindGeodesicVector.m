function [v,X2n,dist] = FindGeodesicVector(q1,q2,X2,lam,flag)

n = size(q1,1);
if (norm(q1-q2) < 0.001)
    v = zeros(size(q1));
    dist = 0;
    q2n = q2;
    X2n = X2;
else
    A = q1*q2';
    [U,S,V] = svd(A);
    if det(A)> 0
        Ot = U*V';
    else
        if n == 2
            Ot = U*([V(:,1) -V(:,2)])';
        else
            Ot = U*([V(:,1) V(:,2) -V(:,3)])';
        end
    end
    X2 = Ot*X2;
    q2 = Ot*q2;
    
    if(flag)
        [gam] = DynamicProgrammingQ(q2/sqrt(InnerProd_Q(q2,q2)),q1/sqrt(InnerProd_Q(q1,q1)),lam,0);
%         [gam] = DynamicProgrammingQ(q2,q1,lam,0);
        gam = (gam-gam(1))/(gam(end)-gam(1));
        X2n = Group_Action_by_Gamma_Coord(X2,gam);
        q2n = curve_to_q(X2n);
    else
        X2n = X2;
        q2n = q2;
    end
    
    N = size(X2n,2);
    dist = real(acos(trapz(sum(q1.*q2n))/(N-1)));
%     dist = InnerProd_Q(q1-q2,q1-q2);
    
    v = q2n - q1;
    v = v - (trapz(sum(v.*q1))/(N-1))*q1;
    nv = sqrt(trapz(sum(v.*v))/(N-1));
    v = dist*v/nv;
end