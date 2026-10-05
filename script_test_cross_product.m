

syms a1 a2 a3 b1 b2 b3

a = [a1;a2;a3];b=[b1;b2;b3];


K = @(a) [0 -a(3) a(2); a(3) 0 -a(1);-a(2) a(1) 0];

KA = K(a)
KB = K(b)

a*b.' - b*a.'

KAB = K(cross(a,b))