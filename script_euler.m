clear

syms phi theta psi r1 r2 r3;
c1 = cos(phi);
s1 = sin(phi);
c2 = cos(theta);
s2 = sin(theta);

c3 = cos(psi);
s3 = sin(psi);

R = [c2*c3 c2*s3 -s2;
    -c1*s3+s1*s2*c3 c1*c3+s1*s2*s3 s1*c2;
    s1*s3 + c1*s2*c3 -s1*c3 + c1*s2*s3 c1*c2;].'
D1 = diff(R,phi)
D2 = diff(R,theta)
D3 = diff(R,psi)

rbar = [r1;r2;r3]
D= [D1*rbar D2*rbar D3*rbar]

Dt = [D1 D2 D3]*[rbar rbar rbar]

simplify(D-Dt)