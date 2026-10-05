clear;

psi = 10;

%clockwise rotation
R = [cosd(psi) sind(psi);-sind(psi) cosd(psi)];

t = [4;5];

p = [6;3];

p1 = [10;0];
p2 = [0;10];

r1 = R*p1*2/3 + t;
r2 = R*p2*2/3 + t;


r = R*p + t;

figure(1);
subplot(1,3,3);

plot([t(1) r1(1)],[t(2) r1(2)],'LineWidth',2,'Color','r');hold on;axis equal;
plot([t(1) r2(1)],[t(2) r2(2)],'LineWidth',2,'Color','r');grid on;

plot([0 p1(1)],[0 p1(2)],'LineWidth',2,'Color','b');
plot([0 p2(1)],[0 p2(2)],'LineWidth',2,'Color','b');

plot([0 t(1)],[0 t(2)],'LineWidth',2,'Color','k');
plot([t(1) r(1)],[t(2) r(2)],'LineWidth',2,'Color','g')
plot([0 r(1)],[0 r(2)],'LineWidth',2,'Color','m')



px = [p(1);0];
py = [0;p(2)];

rx = R*px + t;
ry = R*py + t;

plot([rx(1) r(1)],[rx(2) r(2)],...
    'LineWidth',2,'Color','g','LineStyle','--');
plot([ry(1) r(1)],[ry(2) r(2)],...
    'LineWidth',2,'Color','g','LineStyle','--');
plot([r(1) r(1)],[0 r(2)],...
    'LineWidth',2,'Color','m','LineStyle','--');
plot([r(1) 0],[r(2) r(2)],...
    'LineWidth',2,'Color','m','LineStyle','--');

axis


