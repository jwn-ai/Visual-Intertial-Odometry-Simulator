clear;

syms Px Pv Pxv Qx Qv r0 v0 R  T  Qa




x1 = [0 0 0].';
F = [1 T 0;0 1 T;0 0 1];
P1 = 0*eye(3);
Q = [Qx 0 0;0 Qv 0;0 0 Qa];
H = [1 0 0];
I = eye(3);

%n=1
syms z1 Kx1 Kv1 Ka1

x0 = F*x1

P0 = F*P1*F.' + Q

S = H*P0*H' + R
Sinv = inv(S);
i = z1 - H*x0
%K = P0*H.'*Sinv
K = [Kx1; 0;0]
x1 = x0 + K*i
P1 = (I - K*H)*P0

%n=2;
syms z2 Kx2 Kv2 Ka2

x0 = F*x1

P0 = F*P1*F.' + Q

S = H*P0*H' + R
Sinv = inv(S);
i = z2 - H*x0
%K = P0*H.'*Sinv
K = [Kx2; Kv2; 0]
x1 = x0 + K*i
P1 = (I - K*H)*P0

%n=3;
syms z3 Kx3 Kv3 Ka3

x0 = F*x1

P0 = F*P1*F.' + Q

S = H*P0*H' + R
Sinv = inv(S);
i = z3 - H*x0
%K = P0*H.'*Sinv
K = [Kx3; Kv3;Ka3]
x1 = x0 + K*i
P1 = (I - K*H)*P0





