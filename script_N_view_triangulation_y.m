clear;
%optical axis is along y
J = 3;
B = 1000;
fl = 10;
K = [fl 0 0;0 fl 0;0 0 1];

rn = [B*(0:J-1);0*rand(1,J)*B/2;0*rand(1,J)*B/6];
T = 1000;
phi = 0;
theta = 0;
psi = pi/10; %z
r = [B*(J-1)/2+50;2*B;0];

for i=1:J
    R = utils.e2R(phi,theta,psi); %r = R*p + t (R: C to W)
    RR(:,:,i) = [1 0 0;0 0 1;0 1 0]*R;
    P(:,:,i) = K*[R.' -R.'*rn(:,i)];
end


sigma = 0:0.5:2;

for n=1:length(sigma)
    for t=1:T
        for i=1:J

            p(:,i) = utils.TF(r,RR(:,:,i),rn(:,i));

            y(:,i) = (utils.PIz(p(:,i),fl)) + sigma(n)*randn(2,1);

        end

        rhat = PnP(P,y);

        rhat2 = utils.TriangulateLS(y,P);

        e1(t) = 1/sqrt(3)*norm(rhat - r);
        e2(t) = 1/sqrt(3)*norm(rhat2 - r);

    end

    E1(n) = mean(e1);
    E2(n) = mean(e2);

end

figure,plot(sigma,E1,'LineWidth',2);hold on;grid
plot(sigma,E2,'LineWidth',2);
ylabel('RMSE/B');xlabel('Pixel STD');
title('DLT Comparison')
legend('Homogenous LS','Standard LS')

function rhat = PnP(P,y)
H = [];
N = size(y,2);
for n=1:N
    xn = y(1,n);
    yn = y(2,n);

   Pn = P(:,:,n);
    a1 = Pn(1,:);
    a2 = Pn(2,:);
    a3 = Pn(3,:);
    Hn = [a3*xn - a1; a3*yn-a2];
    H = [H; Hn];

    

end
[U,S,V] = svd(H);
    rhat = V(:,end);
    s = rhat(4);

        rhat = rhat/s;
        rhat = rhat(1:3);
end