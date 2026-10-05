clear;close all;
imu_header;

Target_data = 'target.mat';
%camera optical axis pointing in direction of travel
%attitude 

[r,t,q,RR,w,w2b,v,v2,v2b,a,a2b] = RobotPath(T);
N = length(r);
v0 = gradient(r,T); %central difference to time match derivatives
u = v0./vecnorm(v0,2);



a0 = gradient(v0,T);
%w = cross(r./sqrt(sum(r.^2,1)),v); %particle w wrt origin
%w(isnan(w)) = 0;
x = [r;v0;q;zeros(3,size(r,2));zeros(3,size(r,2))];

%define body frame
for n=1:N
    R = RR(:,:,n);
    u1(:,n) = R(:,1);
    u2(:,n) = R(:,2);
    u3(:,n) = R(:,3);
    ab(:,n) = R.'*a0(:,n);
    wb(:,n) = w(:,n);
    M(:,:,n) = [R.' -R.'*r(:,n)];  %rotation/trans for each view
    Pcam_project(:,:,n) = L*Rbc.'*M(:,:,n); %CAMERA projection matrix with axis swap
end


%
%Simulated IMU inputs
z0 = [ab;wb];

rr(:,1) = r(:,1);rr0(:,1)=r(:,1);
v2(:,1) = v0(:,1);
for n=2:N
    v2(:,n) = v2(:,n-1) + RR(:,:,n)*(ab(:,n))*T;

rr0(:,n) = rr0(:,n-1) + v0(:,n)*T+RR(:,:,n)*ab(:,n)*T^2/2;
rr(:,n) = rr(:,n-1) + v2(:,n)*T;

end
figure,plot(v0(1,:));hold on;plot(v(1,:));plot(v2(1,:));

figure,plot(x(1,:));hold on;plot(r(1,:));plot(rr(1,:));plot(rr0(1,:));
figure,plot(x(2,:));hold on;plot(r(2,:));plot(rr(2,:));
figure,plot(x(3,:));hold on;plot(r(3,:));plot(rr(3,:));
figure,plot(x(4,:));hold on;plot(v2(1,:));plot(v0(1,:))
figure,plot(x(5,:));hold on;plot(v2(2,:));plot(v0(2,:))
figure,plot(x(6,:));hold on;plot(v2(3,:));plot(v0(3,:))

opts.Kmin = 6;
opts.maxLandmarks = 100;
opts.anchorStride = 25;
opts.coverageStride = 5;
opts.candidatesPerAnchor = 30;
opts.depthRange = [5 20];
opts.depthMin = 1.0;
opts.marginPx = 100;

[LM, S,mapInfo] = buildTrajectoryAwareMap( ...
    RR, x(1:3,:), Rbc, fl, imageWidthPx, imageHeightPx, opts);
save(Target_data ,'x','z0','t','T','RR','S','LM','M','Pcam_project','mapInfo');

figure;
plot(mapInfo.coverage);
yline(opts.Kmin,'r--');
ylabel('Visible landmarks');
xlabel('Time step');





% if 0%isfile(Target_data)
%     load('target.mat');
% else
%     LM = GenLandMarks(r,u2,phi); 
%     %Y = Camera(r,[0;1;0],fl,RR,phi,l);
%     S = Camera1(r,[0;1;0],LM,fl,RR,phi);
%     %N0 = 4*750;K0 = size(LM,2);
%    %[S,LM] = Ytest(S,N0,K0,LM);
%    %z0 = z0(:,1:N0);
% end


orange = [1,0.5,0];
figure;
hr = plot3(r(1,:),r(2,:), r(3,:));hold on;
hm = plot3(LM(1,:),LM(2,:),LM(3,:),'mx');

as = .5;
for i=2:1/T:length(r)

    h1 = quiver3(r(1,i),r(2,i),r(3,i),u1(1,i),u1(2,i),u1(3,i),'Color','g',...
        'LineWidth',2,'AutoScaleFactor',as);

    h2 = quiver3(r(1,i),r(2,i),r(3,i),u2(1,i),u2(2,i),u2(3,i),...
        'r','LineWidth',2,'AutoScaleFactor',as);
    h3 = quiver3(r(1,i),r(2,i),r(3,i),u3(1,i),u3(2,i),u3(3,i),...
        'b','LineWidth',2,'AutoScaleFactor',as);
    hv = quiver3(r(1,i),r(2,i),r(3,i),u(1,i),u(2,i),u(3,i),'Color','y',...
        'LineWidth',2,'AutoScaleFactor',as,'LineStyle','--');
end

legend([h1, h2, h3, hv],{'xb','yb','zb','Velocity'})
grid on;axis equal;xlabel('x');ylabel('y');zlabel('z');



%estimate landmark locations via least squares from known robot position
% k0 = 7;
% LM(:,k0)
% rhat2 = utils.TriangulateLS(S(:,1:4,k0),Pcam_project);
% % 
%  Xn(:,1) = x(1:3,1);
% for n=2:N-1
% 
%     xn = testk(x(1:3,n-1),v0(:,n-1),v0(:,n),a0(:,n-1),a0(:,n),[0;0;0],T,eye(3));
%     Xn(:,n) = xn;
% end
% figure;
% plot3(Xn(1,:),Xn(2,:), Xn(3,:));hold on;title('path')
% plot3(r(1,:),r(2,:), r(3,:));hold on;title('path');axis equal;

%--------------------------------------------------------------------





function LM = GenLandMarks(r,u,phi)

N = 1;

Nv = length(u);
f = [2;2;2];
LM = zeros(3,Nv*N);  
k=0;
for n=1:Nv
    rn = r(:,n);
    un = u(:,n);
    i=0;
    while i<N
        p = rn + f.*rand(3,1)-f/2;
        d = un.'*(p-rn);
        theta = acos(d/norm(p-rn));
        if abs(theta)<phi/2 && d>0
            k = k + 1;
            LM(:,k) = p;
            i = i + 1;          
        end

    end
end

end

function [r,t,q,RR,w,w2b,v,v2,v2b,a,a2b] = RobotPath(dt)

s1 = 5;
T1 = 1;
t1 = 0:dt:T1;
n1 = T1/dt+1;

ax=0;ay=0;az=11;
r0 = [0;0;0];
g = 10;
theta = pi/4;
psi = pi/4;
a1 = [ax;ay;az-g];
u1 = [cos(theta)*cos(psi); cos(theta)*sin(psi);sin(theta)];
v1 = s1*u1;


u2 = [0 0 1].';
T2 = 10;
t2 = 0:dt:T2;
theta2 = 0.7*2*pi;
R = 2;
r20 = [R;0;0];
omega = u2*theta2/T2;
s2 = R*omega;
n2 = T2/dt+1;




r1 = r0 + v1*t1+.5*a1*t1.^2;

[r2,qq2,RR] = qRotate(u2,t2,norm(omega),r20);


v2b = cross(repmat(omega,1,n2),r2);
w2b = cross(r2./vecnorm(r2).^2,v2b);
a2b = norm(omega)^2*R*[-1;0;0];
r2 = r2-r2(:,1) + r1(:,end);

v2 = cross(repmat(omega,1,n2),r2);
a2 = cross(repmat(omega,1,n2),v2);

v = [v1*ones(1,n1-1)+a1*t1(1:end-1) v2];
r = [r1(:,1:end-1) r2];
t = 0:dt:(T1+T2);
q = [repmat([1;0;0;0],1,n1-1) qq2];
RR = cat(3,repmat(eye(3),1,1,n1-1), RR);
w = [zeros(3,n1-1) omega*ones(1,n2)];

a = [a1*ones(1,n1-1) a2];

end

function Y = CameraSim(r,RR,Rbc,Kmin,imageWidthPx,imageHeightPx,L,Fs,ns,z0,sigma_z,marginPx)
%ensure keyframes have at least Kmin features visible
N = size(r,2);
Y = zeros(2,N,K);
I = zeros(Kmin,ceil(N/ns));
while c<Kmin
    c = 0;
    %generate test sample point
    y = [imageWidthPx/2 0; 0 imageHeightPx/2]*(2*rand(2,1)-1) * 0.7;
    z = z0 + sigma_z*rand;  % choose useful scene depth [5,20] m
    pc = utils.invPz(y,L,z);
    for n=1:ns:N
        rn = r(:,n);
        R = RR(:,:,n)*Rbc;


        % Sample safely inside the sensor, in pixels


        pw = R*pc + rn;
        I1 = pc(3) > depth_min ;
        I2 =  abs(y(1)) <  imageWidthPx/2-marginPx;
        I3 = bs(y(2)) < imageHeightPx/2-marginPx;
        if (I1  && I2 && I3)

            I(k, n) = 1;
            c = c  + 1;
        end

    end
    k = k + 1;



end
end

function Y = Camera(r,ez,f,RR,phi,l)
%generate LM in fov of all vehicle positions
%1)get ijth local u,v
%2)pick, Z project to local XYZ
%3)transform to nth view
%4) check FOV
N = size(r,2);
M = 10;
a = l/M;
k = 0;
N0 = 25;
for n=1:N0
    n
    rn = r(:,n);
    R = RR(:,:,n);
    for i=1:M
        u = a*(i - M/2);
        for j=1:M
            v = a*(j - M/2);
            y = [u;v];
            z = 10;
            pk = utils.invPIz(y,f,z); %to XYZ
            rk = utils.FF(pk,R,rn);%to global
            cnt = 0;
            for m = 1:N0
                rm = r(:,m);
                Rm = RR(:,:,m);
                pm = utils.TF(rk,Rm,rm);
                if utils.CheckFOV(pm,0,ez,phi)
                    cnt = cnt + 1;
                end
            end
            if cnt==N0
                k = k+1;
                Y(:,n,k) = y;
            end
        end
    end
end
end
   

function [Y1,LM1] = Ytest(Y,N0,K0,LM)

[~,N,K]=size(Y);
I = zeros(K0,N);
    for k = 1:K0

        for n=1:N
            disp(n)
            y = Y(:,n,k);
            if sum(abs(y))>0
          
                I(k,n) = 1;
            end
        end
    end

I   
j=0;
for k=1:K0
    if sum(I(k,1:N0))==N0
        j = j + 1
        for n=1:N0
            Y1(:,n,j) = Y(:,n,k);
            LM1(:,j) = LM(:,k);
        end
    end
end

end

