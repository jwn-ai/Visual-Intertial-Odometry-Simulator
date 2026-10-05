function [r,qq,RR] = qRotate(u,t,omega,r0)

r(:,1) = r0;
qq(:,1) = [1;0;0;0];
RR(:,:,1) = eye(3);
for n=2:length(t)
    mu = omega*t(n);
    q = [cos(mu/2);u*sin(mu/2)];
    qv = q(2:4);
    I = eye(3);

    q1 = q(2);q2 = q(3);q3=q(4);q0 = q(1);
    K = [0 -q3 q2; q3 0 -q1;-q2 q1 0];

    a = q0^2 - (q1^2 + q2^2 + q3^2);
    %cumulative orientation change
    %R is the rotation operating on body frame vectors
    %it is the to matrix R

    R = a*I + 2*(qv*qv.') + 2*q0*K; 

    r(:,n) = R*r0;
    qq(:,n) = q;
    RR(:,:,n) = R;
end

end