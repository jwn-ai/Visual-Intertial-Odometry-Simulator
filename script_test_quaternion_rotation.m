clear;







u = [0 0 -1].';

t1 = 10;
omega = 2*pi/t1;
t = 0:.1:t1;
r0 = [10 0 10].';

r(:,1) = r0;

for n=2:length(t)
    mu = omega*t(n);
q = [cos(mu/2);u*sin(mu/2)];
v = q(2:4);
I = eye(3);

q1 = q(2);q2 = q(3);q3=q(4);q0 = q(1);
K = [0 -q3 q2; q3 0 -q1;-q2 q1 0];

a = q0^2 - (q1^2 + q2^2 + q3^2);

R0 = [q0^2 + q1^2-q2^2-q3^2 2*(q1*q2-q0*q3) 2*(q1*q3+q0*q2)
    2*(q1*q2+q0*q3) q0^2 - q1^2+q2^2-q3^2 2*(q2*q3-q0*q1)
    2*(q1*q3-q0*q2) 2*(q2*q3+q0*q1) q0^2 - q1^2-q2^2+q3^2];


%a = (2*q0^2-1)
R = a*I + 2*(v*v.') + 2*q0*K;

r(:,n) = R*r0;

end


figure;
plot3(r(1,:),r(2,:),r(3,:),'.');


mu = pi/4;
q = [cos(mu/2);u*sin(mu/2)];
v = q(2:4);
I = eye(3);

q1 = q(2);q2 = q(3);q3=q(4);q0 = q(1);
K = [0 -q3 q2; q3 0 -q1;-q2 q1 0];

a = q0^2 - (q1^2 + q2^2 + q3^2);
R = a*I + 2*(v*v.') + 2*q0*K;
