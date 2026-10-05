

clear;

%v = v + R(q)(ahat)
syms q0 q1 q2 q3 r0 r1 r2  w;
v = [q1;q2;q3];

q = [q0 v.'].';

I = eye(3);

K = [0 -q3 q2; q3 0 -q1;-q2 q1 0];

d = q0^2 - (q1^2 + q2^2 + q3^2);
%a = (2*q0^2-1)
R = d*I + 2*(v*v.') + 2*q0*K


r = [r0;r1;r2];


%for Rtranspose
Rt = R.';
D0 = diff(Rt,q0)
D1 = diff(Rt,q1)
D2 = diff(Rt,q2)
D3 = diff(Rt,q3)

Dt = [D0*r D1*r D2*r D3*r]

Dt2 = 2*[q0*r0 - q2*r2 + q3*r1, q1*r0 + q2*r1 + q3*r2, q1*r1 - q0*r2 - q2*r0, q0*r1 + q1*r2 - q3*r0;
q0*r1 + q1*r2 - q3*r0, q0*r2 - q1*r1 + q2*r0, q1*r0 + q2*r1 + q3*r2, q2*r2 - q0*r0 - q3*r1;
q0*r2 - q1*r1 + q2*r0, q3*r0 - q1*r2 - q0*r1, q0*r0 - q2*r2 + q3*r1, q1*r0 + q2*r1 + q3*r2]
 
Dt - Dt2

%from matlab

Dx = 2*[q0*r0 - q2*r2 + q3*r1,  q1*r0 + q2*r1 + q3*r2, q1*r1 - q0*r2 - q2*r0, q0*r1 + q1*r2 - q3*r0;
        q0*r1 + q1*r2 - q3*r0,  q0*r2 - q1*r1 + q2*r0, q1*r0 + q2*r1 + q3*r2, q2*r2 - q0*r0 - q3*r1;
        q0*r2 - q1*r1 + q2*r0,  q3*r0 - q1*r2 - q0*r1, q0*r0 - q2*r2 + q3*r1, q1*r0 + q2*r1 + q3*r2];

Dx-Dt2