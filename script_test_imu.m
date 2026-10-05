
clear;close all;

Tr = 5;  %circle rotation time
u = [0 0 1];
wr = 2*pi/Tr;  %angular speed
dt = 1/100;
t = 0:dt:Tr;
z0 = 10;
x0 = [1 0 0].';
X(:,1) = [0 0 0].';
l = 1/10;
I = eye(3);
R2(:,:,1) = eye(3);
w = wr*u;
W = -[0 w(3) -w(1);-w(3) 0 w(2);w(1) -w(2) 0];

for n=2:length(t)
    tn = t(n);
mu = wr*tn;  %angle vs time
R = [cos(mu) -sin(mu) 0; sin(mu) cos(mu) 0; 0 0 1]; %rotation around z
Rs(:,:,n) = R;
R2(:,:,n) = (I+dt*sinc(wr*dt)*W + .5*dt^2*(sinc(wr*dt/2))^2*W^2)*R2(:,:,n-1);
x = R*x0;
X(:,n) = x;
plot3([X(1,n)], [X(2,n)],[X(3,n)],'.');hold on;axis equal;grid;


col = ['rgb'];
if mod(n,10)==0

    for i=1:3
        u = R(:,i)*l;
        quiver3(x(1),x(2),x(3),u(1),u(2),u(3),col(i));hold on;
    end
end




end


figure

    for i=1:3
         u = R(:,i)*l;
        quiver3(x(1),x(2),x(3),u(1),u(2),u(3),col(i));hold on;
    end