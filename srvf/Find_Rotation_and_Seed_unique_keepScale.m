function [q2best,Rbest] = Find_Rotation_and_Seed_unique_keepScale(q1,q2,reparamFlag)
% This function returns a locally optimal rotation and seed point for shape
% q2 w.r.t. q1

[n,T] = size(q1);

scl = 10;
minE = 1000;
for ctr = 0:floor((T)/scl)
    q2n = ShiftF(q2,scl*ctr);
    [q2n,R] = Find_Best_Rotation(q1,q2n);
    
    if(reparamFlag)
        
        if norm(q1-q2n,'fro') > 0.0001
            gam = DynamicProgrammingQ(q1,q2n,0,0);
            gamI = invertGamma(gam);
            gamI = (gamI-gamI(1))/(gamI(end)-gamI(1));
            p2n = q_to_curve(q2n);
            p2new = Group_Action_by_Gamma_Coord(p2n,gamI);
            q2new = curve_to_q_keepScale(p2new);
            q2new = ProjectC_keepScale(q2new);
        else
            q2new = q2n;
        end
        
    else
        q2new  = q2n;
    end
    %Ec = acos(InnerProd_Q(q1,q2new));
    %Ec = InnerProd_Q(q1-q2new,q1-q2new);
    Ec = norm(q1-q2new,'fro');
    
    if Ec < minE
        Rbest=R;
%         Rbest=[1,0;0,1];
        q2best  = q2new;
        minE = Ec;
    end
end

return;