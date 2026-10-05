clear;

syms q0 q1 q2 q3 r0 r1 r2 w;
v = [q1;q2;q3];

q = [q0 v.'].';

I = eye(3);

K = [0 -q3 q2; q3 0 -q1;-q2 q1 0];

a = q0^2 - (q1^2 + q2^2 + q3^2);
%a = (2*q0^2-1)
R = a*I + 2*(v*v.') + 2*q0*K

R2 = I +2*q0*K+2*K^2

r = [r0;r1;r2];
%for R
D0 = diff(R,q0)
D1 = diff(R,q1)
D2 = diff(R,q2)
D3 = diff(R,q3)

D = [D0*r D1*r D2*r D3*r]

s2 = [q0;-q3;q2];
s3 = [q3;q0;-q1];
s4 = [-q2;q1;q0];
s1 = -[q1;q2;q3];

S = [s1.';s2.';s3.';s4.']
A = S*r

a1 = s1.'*r;
a2 = s2.'*r;

a3 = s3.'*r;

a4 = s4.'*r;

D2 = [a2 -a1 a4 -a3;a3 -a4 -a1 a2;a4 a3 -a2 -a1]

D3 = [A(1) A(4) -A(3) A(2);A(2) A(3) A(4) -A(1);A(3) -A(2) A(1) A(4)]


%for Rt
Rt = R.';
D0 = diff(Rt,q0)
D1 = diff(Rt,q1)
D2 = diff(Rt,q2)
D3 = diff(Rt,q3)

Dt = [D0*r D1*r D2*r D3*r]



%another identity for D
w = cross(v,r);
K = @(a) [0 -a(3) a(2); a(3) 0 -a(1);-a(2) a(1) 0];

D4 = -[w-q0*r (-q0*K(r)+K(w)) - I*(r.'*v)]


e1 = I(:,1);e2 = I(:,2);e3=I(:,3);
K1 = K(e1);K2=K(e2);K3=K(e3);

2*(q0*K1 + (K1*K(v) + K(v)*K1))

D5 = (q0*I+K(v))*K1*R


I4 = eye(4);
I3 = eye(3);


syms q0 q1 q2 q3 w1 w2 w3;
v = [q1;q2;q3]
q = [q0 v.'].'



Omega = @(v) [0 -v(3) v(2);v(3) 0 -v(1);-v(2) v(1) 0 ];

ML = @(q) q(1)*I4 + [0 -v.';v Omega(v)];
MR = @(q) q(1)*I4 + [0 -v.';v -Omega(v)];

M2 = @(q) [-v.';q(1)*I3 + Omega(v)];
M3 = @(q) [q M2(q)];


w = [w1;w2;w3];

ML(q)*[0;w]

M2(q)

M3(q)


