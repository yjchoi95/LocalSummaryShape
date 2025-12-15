function d = GeodesicElasticClosed(p1,p2,dt)

% input p1 and p2 as 2xn matrices
% to turn off figures set figs=0
Disp_registration_between_curves=1;
figs=1;

stp = 6;

if figs
figure(1); clf; hold on;
plot(p1(1,:),p1(2,:),'r','LineWidth',2);
plot(p2(1,:),p2(2,:),'b','LineWidth',2);
axis equal;
axis xy off;
end

p1 = ReSampleCurve(p1,100);
p2 = ReSampleCurve(p2,100);

q1 = curve_to_q(p1);
q2 = curve_to_q(p2);

tic
[q2n,~] = Find_Rotation_and_Seed_unique(q1,q2,1);
q2n = q2n/sqrt(InnerProd_Q(q2n,q2n));
q2n = ProjectC(q2n);
p2n = q_to_curve(q2n);
toc

% d = acos(InnerProd_Q(q1,q2n));
d = sqrt(InnerProd_Q(q1-q2n,q1-q2n));

if figs
alpha = geodesic_sphere_Full(q1,q2n,stp);
for i = 2:(size(alpha,3)-1)
    alpha(:,:,i) = alpha(:,:,i)./sqrt(InnerProd_Q(alpha(:,:,i),alpha(:,:,i)));
end
% alpha1 = geodesic_sphere_Full_noproj(q1,q2n,stp);
Path_Plot(alpha,p2n,[73,6],dt);
axis xy; view([0 90]);
% Path_Plot(alpha1,11,'b',[73,6]);
% axis xy;
X1 = q_to_curve(q1);
X2 = q_to_curve(q2n);
figure;hold on;
plot([X1(1,:) X1(1,1)],[X1(2,:) X1(2,1)],'-','LineWidth',2);
plot([X1(1,:) X1(1,1)],[X1(2,:) X1(2,1)],'*','LineWidth',2);
plot(X1(1,1),X1(2,1),'*r','LineWidth',10);
plot(X1(1,25),X1(2,25),'*g','LineWidth',10);
plot(X1(1,50),X1(2,50),'*k','LineWidth',10);
plot(X1(1,75),X1(2,75),'*c','LineWidth',10);
axis equal off; view([1 90])
figure(3),clf,hold on;
plot([X2(1,:) X2(1,1)],[X2(2,:) X2(2,1)],'-','LineWidth',2);
plot([X2(1,:) X2(1,1)],[X2(2,:) X2(2,1)],'*','LineWidth',2);
plot(X2(1,1),X2(2,1),'*r','LineWidth',10);
plot(X2(1,25),X2(2,25),'*g','LineWidth',10);
plot(X2(1,50),X2(2,50),'*k','LineWidth',10);
plot(X2(1,75),X2(2,75),'*c','LineWidth',10);
axis equal off; view([1 90])
end



if figs && Disp_registration_between_curves
% Displaying the correspondence
    %X2n=alpha(:,:,end);
    %X1=alpha(:,:,1);
    X2n=X2;
    figure; clf;
    z = plot(X1(1,:), X1(2,:),'r-+');
    set(z,'LineWidth',2);
    axis off;
    hold on;
    z = plot(X2n(1,:), 0.25+X2n(2,:),'b-+');
    set(z,'LineWidth',2);
    N = size(X1,2);
    for j=1:N/15
        i = j*15;
        plot([X1(1,i) X2n(1,i)],[X1(2,i) 0.25+X2n(2,i)], 'k','Linewidth',1.5);
    end
end