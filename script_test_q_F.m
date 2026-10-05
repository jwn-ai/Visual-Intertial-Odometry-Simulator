clear;

%target trajectory gives global a and w.
%with F and a,v,w should be able to duplicate track
%a and w are converted to body frame
%use these to regenerate target track using F matrix


load target.mat;



I = eye(3);
gamma = 0;
alpha = 0;
N = size(x,2);



x_n(:,1) = x(:,1);
for n=2:N-1
    % Rn = R(:,:,n);
    % ab_n = R_n*a(:,n);
    % wb_n = R_n*w(:,n);
    xn = x(:,n);
    zn = z(:,n);
    an = zn(1:3);
    wn = zn(4:6);
    Rn = RR(:,:,n);
    ab_n = Rn.'*an; %second motion is already in body frame
    wb_n = wn; %already in body frame
    zb_n = [ab_n;wb_n];
    Fn = F(xn,zb_n,T);

    x_n = f(xn,zb_n,T);

    x2(:,n) = x_n;
    x3(:,n) = Fn(xn,zb_n,T);

end

figure,plot3(x2(1,:),x2(2,:), x2(3,:));hold on;title('path')
grid on;axis equal;

function Fx = F(xn,zn,T) 

q = xn(7:10);
a = zn(1:3);
w = zn(4:6);
I = eye(3);
O = zeros(3);
O43 = zeros(4,3);
O34 = zeros(3,4);



Omega = @(v) [0 -v(3) v(2);v(3) 0 -v(1);-v(2) v(1) 0 ];
Gamma = @(q0,v) q0*eye(4) + [0 -v.';v Omega(v)];
Xi = @(q) [-q(2:4).';q(1)*eye(3) + Omega(q(2:4))];

D = @(q,r) [q(1)*r(1) + q(3)*r(3) - q(4)*r(2) q(2)*r(1) + q(3)*r(2) + q(4)*r(3) q(1)*r(3) + q(2)*r(2) - q(3)*r(1) q(2)*r(3) - q(1)*r(2) - q(4)*r(1);
q(1)*r(2) - q(2)*r(3) + q(4)*r(1) q(3)*r(1) - q(2)*r(2) - q(1)*r(3) q(2)*r(1) + q(3)*r(2) + q(4)*r(3), q(1)*r(1) + q(3)*r(3) - q(4)*r(2);
q(1)*r(3) + q(2)*r(2) - q(3)*r(1) q(1)*r(2) - q(2)*r(3) + q(4)*r(1), q(4)*r(2) - q(3)*r(3) - q(1)*r(1) q(2)*r(1) + q(3)*r(2) + q(4)*r(3)];
 
R = @(q) (2*q(1)^2-1)*I + 2*(q(2:4)*q(2:4).') + 2*q(1)*Omega(q(2:4));

F1 = [I I*T O34 O O];
F2 = [O I*T D(q,a)*T -R(q)*T O];
F3 = [O43 O43 .5*T*Gamma(0,w) O43 -.5*T*Xi(q)];
F4 = [O O O34 I O];
F5 = [O O O34 O I];
Fx = [F1;F2;F3; F4;F5];
end

function x = f(xn,zn,T)


r = xn(1:3);
v = xn(4:6);
a = zn(1:3);
w = zn(4:6);
q = xn(7:10);
alpha = xn(11:13);
gamma = xn(14:16);


I = eye(3);
Omega = @(v) [0 -v(3) v(2);v(3) 0 -v(1);-v(2) v(1) 0 ];
Gamma = @(q0,v) q0*eye(4) + [0 -v.';v Omega(v)];

R = @(q) (2*q(1)^2-1)*I + 2*(q(2:4)*q(2:4).') + 2*q(1)*Omega(q(2:4));

r = r + v*T +.5*a*T^2;
v = v + R(q)*(a-alpha)*T;
q = .5*Gamma(0,w-gamma)*q;



x = [r;v;q;alpha;gamma];

end
