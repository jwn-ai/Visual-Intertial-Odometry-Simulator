%total state EKF, x = [s m].  s = IMU, m = LM/features
%only observed features used and constant number K per time step
%m and s are separately computed
%C = Pbar (predicted error covariance
% NOTE: comment dimensions are really times 3, e.g. 16 = 3*16
%P = Phat (filtered error cov)
%xhat = xfiltered
%mu = xpredicted

%IMU used to propagate states
%visual measurements used to estimate IMU biases and correct

clear;
close all;
imu_header;
load target.mat;

[~,N,K] = size(S);
Kmin = 6;
Ns = 16;  %16 states in quaternion attitude formulation
Ns1 = Ns-1;
Nsm = Ns+3*K; %nominal dimension [r v q ba bg r1 ... rK]
Nc = Nsm-1;  %correction dimension (angles, not quaternions)

%feature noise parameters
V0 = 0.6;  %pixel errror in 
V = V0^2*eye(2);
Pm0 = (1)^2;
ny = V0*randn(size(S));
Y = S + ny;
P0 = [1e-2 1e-5 1e-4 1e-4 1e-7 Pm0];
Q0 = blkdiag([1e-3 1e-2 2e-4 1e-5 1e-4]);

P = initP(P0,K);

%noisy IMU inputs
nz = chol(Sigma_z,'lower') * randn(6,N);
b = diag([alpha; gamma])*ones(6,N);
z = z0 + b + nz;


%initialize x = [s;m], s=body,m=landmarks

Ni = 70;
for k=1:K
    m(1+3*(k-1):3*k,1) = LM(:,k);
end
mhat = initLM(Y,mapInfo,Pcam_project,Ni); %using truth camera poses for simplicity
sum(abs(m-mhat))
xshat = x(:,1);

%allocate memory for saved data
Mu(:,1) = xshat;
Mu1(:,1) = xshat;
mu0 = xshat;
mu = zeros(Ns,1);
xhat = [xshat;mhat];
PPhat = zeros(Nc,Nc,N); %save data
Xhat = zeros(Nsm,N);
Xhat(:,1) = xhat;
e = zeros(Nsm,N);
d = zeros(N,K);
Ytilde = zeros(2,N,K);

for n=2:N
    %% predict states based on IMU data and motion model
    xshat = xhat(1:Ns);
    rn = xshat(1:3);  %estimated position at t=n-1
    Rn0 = RR(:,:,n); %true attitude 
    qn = xshat(7:10);
    Rn = utils.q2R(qn);
    alpha_hat = xshat(11:13);
    gamma_hat = xshat(14:16);
    mhat = xhat(Ns+1:end);

    zn = (z(:,n)-[alpha_hat;gamma_hat]);  %noisy biased subtracted IMU output
    wn = zn(4:6);
    mu = f(xshat,zn,T); %predicted state

    [Ps,Psm,Pm] = getP(P);
    Fs = F(xshat,zn,T);
    Qn = Q(Sigma_z,Rn,T);
    Cs = Fs*Ps*Fs.' + Qn; %16x16
    Csm = Fs*Psm; %16x3K
    P = [Cs Csm;Csm.' Pm];    %15x15

    %test kinematics
   
    mu0 = f(mu0,  z0(:,n),T); %can show what happens when bias not accounted
    Mu(:,n) = mu0;


    %% sequential measurement update/correction
    xhatk = [mu;xhat(Ns+1:end)];
    Kn = sum(mapInfo.visibility(:,n));
    visibleIds = find(mapInfo.visibility(:,n));
    for id=1:min(Kmin,Kn)
        k = visibleIds(id);

        yk = Y(:,n,k);  %kth sensor out put in frame n
        yk0 = S(:,n,k); %truth
        rk = xhatk(Ns + 1+3*(k-1):Ns + 3*k);  %estimated LM world position
        Rn = utils.q2R(xhatk(7:10));

        %Predict Measurement
        ykhat = h(rk,Rn*Rbc,xhatk(1:3),fl);
        ykhat0 = h(m(1+3*(k-1):3*k),Rn0*Rbc,x(1:3,n),fl);
        % pkdot = Omega(wn)*pk + vn;

        %[Hs,Hm] = testJh(xsk,ykhat,rk,fl,Rbc); %numerical J

        [Hs, Hm] = J_error(xhatk, rk, Rn, fl, Rbc); %analytical J

        H = zeros(2, Nsm-1);
        H(:,1:15) = Hs;
        H(:,16 + 3*(k-1) : 18 + 3*(k-1)) = Hm;

        Sk = H*P*H.' + V;
        G = P*H.'/Sk;

        ytildek = yk - ykhat;
        ytildek0 = yk0 - ykhat0;
        d(n,k) = ytildek.'*(Sk\ytildek);

        xtildek =G*ytildek; %r v dtheta (3x1,7:9) ba bg
        xhatk(1:3)   = xhatk(1:3)   + xtildek(1:3);
        xhatk(4:6)   = xhatk(4:6)   + xtildek(4:6);
        xhatk(7:10)  = utils.correctq(xtildek(7:9), xhatk(7:10));
        xhatk(11:13) = xhatk(11:13) + xtildek(10:12);
        xhatk(14:16) = xhatk(14:16) + xtildek(13:15);
        xhatk(Ns+1:end) = xhatk(Ns+1:end) + xtildek(Ns:end);

        %update vehicle covariance and kth block,Joseph form
        A = eye(size(P)) - G*H;
        P = A*P*A' + G*V*G';
        Greset = eye(size(P));
        Greset(7:9,7:9) = eye(3) - 0.5 * utils.Omega(xtildek(7:9));

        P = Greset * P * Greset.';
        P = (P + P')/2;

        if ~isreal(P) || ~isreal(xhatk)
            fprintf('WARNING: Imaginary numbers leaked into the filter!\n');
        end
        Ytilde(:,n,k) = ytildek;
        

    end %k

    %update state after landmarks observation complete
    xhat = xhatk;
    
    Xhat(:,n) = xhat;

    PPhat(:,:,n) = P;
    if rem(n,500)==0
        ep = sum(abs(mu(1:3)-xhatk(1:3)));
        disp(['predicted error = ',num2str(ep) '.']);
        disp([num2str(cond(P)),', n = ' num2str(n)]);

        e(:,n) = abs(xhat-[x(1:Ns,n);m]);
        disp(['99% NIS Chi-Square = ',num2str(sum(d(:)<9.21)/numel(d(:)))])

        [~, cholFlag] = chol(P);
        lambda = eig(P);

        fprintf('chol flag: %d, min eig: %.3e, max eig: %.3e\n', ...
            cholFlag, min(lambda), max(lambda));

        fprintf('cond(S): %.3e\n', cond(Sk));
    end

end %n

%% DATA ANALYSIS
%plot covariance
figure;
image(5000*sign(P) .* (abs(P).^.7));
title("Error Covariance Matrix")

%generate uncertainty ellipsoids for states
for i=1:N

[v(:,:,i), lambda] = eigs(PPhat(1:3,1:3,i));
l(:,i) = diag(lambda);
end

%generate uncertainty for landmarks
v_m=zeros(3,3,N,K);
l_m =zeros(3,N,K);
for i=1:N
    for j=1:K
        [v_m(:,:,i,j), lambda_m] = eigs(PPhat(Ns1+1 +(j-1)*3:Ns1+3*j,Ns1+1+(j-1)*3:Ns1+3*j,i));
        l_m(:,i,j) = diag(lambda_m);
    end
end

%extract attitude error
for i=1:N
qHat  = Xhat(7:10,i);
qTrue = x(7:10,i);
dtheta(:,i) = 180/pi*utils.rightAttitudeError(qHat, qTrue);
end

figure;
for j=1:3

    subplot(1,3,j)
    plot(dtheta(j,:));hold on;
    ylabel(['Angle coordinate # ', num2str(j) ' [deg]'])
    xlabel('time');
    title('Attitude Error')
end

%plot error
str1 = {'position','velocity','accel bias', 'gyro bias'};
E = x - Xhat(1:Ns,:);
cnt = 0;
for i=[0 3 10 13]
    cnt = cnt+1;
     figure;
    for j=1:3
       
        subplot(1,3,j)
        plot(E(i+j,1:n));hold on;
        ylabel(['Error coordinate # ', num2str(j)])
        xlabel('time');
        title(str1{cnt})
    end
end

Er = E(1:3,:);
disp(['Cumulative position error = ' num2str(sum(abs(Er(:)))) ' m.'])



figure;
for i=1:3
subplot(1,3,i);
plot(Xhat(i,1:n));hold on;
plot(Mu(i,1:n));
plot(x(i,1:n));
legend('estimated','predicted','true');
ylabel(['position # ', num2str(i)])
xlabel('time');
title('Position')
end
figure;
for i=1:3
subplot(1,3,i);
plot(Xhat(i+3,1:n));hold on;
plot(Mu(i+3,1:n));
plot(x(i+3,1:n));
legend('estimated','predicted','true');
ylabel(['velocity # ', num2str(i)])
xlabel('time');
title('Velocity')
end

figure;
plot3(Xhat(1,1:n),Xhat(2,1:n),Xhat(3,1:n),'r');
hold on;grid on;axis equal;
%plot uncertainty
for i=[1:250:N N]
    [X,Y,Z] = ellipsoid(0, 0,0,...
        sqrt(l(1,i)),sqrt(l(2,i)),sqrt(l(3,i)));
    pts = [X(:), Y(:), Z(:)]';
    pts_rotated = v(:,:,i).' * pts + Xhat(1:3,i); % R is your arbitrary 3x3 rotation matrix

    % 3. Put them back into standard grid format
    X_rot = reshape(pts_rotated(1,:), size(X));
    Y_rot = reshape(pts_rotated(2,:), size(Y));
    Z_rot = reshape(pts_rotated(3,:), size(Z));

    surf(X_rot, Y_rot, Z_rot);

    if i==N
        for k=1:K
            [X,Y,Z] = ellipsoid(0, 0,0,...
                sqrt(l_m(1,i,k)),sqrt(l_m(2,i,k)),sqrt(l_m(3,i,k)));
            pts = [X(:), Y(:), Z(:)]';
            pts_rotated = v_m(:,:,i,k) * pts + Xhat(Ns+1+(k-1)*3:Ns+3*k,i); % R is your arbitrary 3x3 rotation matrix

            % 3. Put them back into standard grid format
            X_rot = reshape(pts_rotated(1,:), size(X));
            Y_rot = reshape(pts_rotated(2,:), size(Y));
            Z_rot = reshape(pts_rotated(3,:), size(Z));
            surf(X_rot, Y_rot, Z_rot);

        end
    end
end
plot3(Mu(1,:),Mu(2,:),Mu(3,:),'g');
plot3(x(1,1:n),x(2,1:n),x(3,1:n),'b');axis equal;
plot3(LM(1,:),LM(2,:),LM(3,:),'co')
plot3(mhat(1:3:end),mhat(2:3:end),mhat(3:3:end),'mx')
legend('estimated','truth')

figure;
plot3(Xhat(1,1:n),Xhat(2,1:n),Xhat(3,1:n),'r');
hold on;grid on;axis equal;
%plot uncertainty
for i=1:250:N
    [X,Y,Z] = ellipsoid(0, 0,0,...
        sqrt(l(1,i)),sqrt(l(2,i)),sqrt(l(3,i)));
    pts = [X(:), Y(:), Z(:)]';
    pts_rotated = v(:,:,i) * pts + Xhat(1:3,i); % R is your arbitrary 3x3 rotation matrix

    % 3. Put them back into standard grid format
    X_rot = reshape(pts_rotated(1,:), size(X));
    Y_rot = reshape(pts_rotated(2,:), size(Y));
    Z_rot = reshape(pts_rotated(3,:), size(Z));

    surf(X_rot, Y_rot, Z_rot);
end
plot3(Mu(1,:),Mu(2,:),Mu(3,:),'g');
plot3(x(1,1:n),x(2,1:n),x(3,1:n),'b');axis equal;
plot3(LM(1,:),LM(2,:),LM(3,:),'co')
plot3(mhat(1:3:end),mhat(2:3:end),mhat(3:3:end),'mx')
legend('estimated','predicted','truth','LM truth','LM estimated')
title('VIO-SLAM EKF Track')
xlabel('X_{world} [m]');ylabel('Y_{world}, [m]');zlabel('Z_{world}, [m]')
LM_label = string(1:K);
text(LM(1,:)+.3,LM(2,:),LM(3,:)+.3,LM_label)

figure
leg_str = '';
cmap = hsv(K);
shuffled_indices = randperm(K);
cmap = cmap(shuffled_indices, :);

for k=1:K
    plot(d(:,k),'Color',cmap(k,:));hold on; grid on;
    leg_str{k} = ['LM ' num2str(k)];
end
legend(leg_str);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Fx = F(xn,zn,T)

%x = [r v mu alpha gamma], [range,velocity,angle error, ba error, bg error]

q = xn(7:10);
a = zn(1:3);
w = zn(4:6);
I = eye(3);
O = zeros(3);


R = utils.q2R(q);

Fr = [ I, I*T, ...
    -0.5*R*utils.Omega(a)*T^2, ...
    -0.5*R*T^2, ...
    zeros(3) ];

Fvmu = -R*utils.Omega(a)*T;

theta = T*norm(w);
Omega = utils.Omega(w*T);
Fmug = -T*(I - ((1-cos(theta))/theta^2)*Omega + ((theta-sin(theta))/theta^3)*Omega^2);
%Fmug = -T*I;
Fv = [O I Fvmu -R*T O];

%Fmumu = eye(3)-utils.Omega(w)*T; %=exp(-Omega(omega*T) = RotationMatrix(omega*T)^T, is a deltaR
Fmumu = utils.v2R(w*T).';
Fmu = [O O Fmumu O Fmug];
Fa = [O O O I O];
Fg = [O O O O I];
Fx = [Fr;Fv;Fmu; Fa;Fg];

end

function x = f(xn,zn,T)


r = xn(1:3);
v = xn(4:6);
a = zn(1:3);
w = zn(4:6);
q = xn(7:10);
alpha = xn(11:13);
gamma = xn(14:16);


Rn = utils.q2R(q);
%propagate q exact 0th order integration
q = utils.Gamma(utils.qwT(w,T))*q; 
%q2 =utils.Gamma(utils.qwT(w,T))*q;

% Rn2 = utils.q2R(q2);
% 
% Rnhat = (Rn2+Rn)/2; %midpoint
r = r + v*T+Rn*a*T^2/2;
v = v + Rn*a*T;

x = [r;v;q;alpha;gamma];

end


function H = Hc(p,f)

z = p(3);x=p(1);y=p(2);
H = zeros(2,3);
H(1,1) = f/z;
H(2,2) = f/z;
H(1,3) = -f*x/z^2;
H(2,3) = -f*y/z^2;
end

function hx = h(rk,R,r,fl)
    pk = utils.TF(rk,R,r);
    hx = utils.PIz(pk,fl);
end




function [Ps,Psm,Pm] = getP(P)
Ps = P(1:15,1:15);
Psm = P(1:15,16:end);
Pm = P(16:end,16:end);
end



function P = initP(Ps0,K)
I3 = eye(3);
Ps = blkdiag(Ps0(1)*I3,Ps0(2)*I3,Ps0(3)*I3,Ps0(4)*I3,Ps0(5)*I3);
Ns = size(Ps,1);
Pm = Ps0(6)*eye(3*K);
Psm = 0*Ps0(6)*eye(Ns,3*K);
P = [Ps Psm; Psm.' Pm];

end





function [Hs,Hm] = J_error(x, rk, R, fl, Rbc)

    r = x(1:3);

    d  = rk - r;
    pb = R.' * d;         % landmark in body coordinatesm, mk
    pc = Rbc.' * pb;      % landmark in camera coordinates, pk

    Jpi = Hc(pc, fl);     % must use CAMERA point

    % Error-state ordering:
    % [dr(1:3), dv(4:6), dtheta(7:9), dba(10:12), dbg(13:15)]
    Hs = zeros(2,15);

    % Vehicle position
    Hs(:,1:3) = -Jpi * Rbc.' * R.';

    % Right-multiplicative attitude injection:
    % R_new = R * (I + Omega(dtheta))
    Hs(:,7:9) = Jpi * Rbc.' * utils.Omega(pb);

    % Landmark-world-coordinate derivative
    Hm = Jpi * Rbc.' * R.';
end

function Qk = Q(Sigma_z,Rn,T)

L = zeros(15,6);

L(1:3,1:3) = -0.5 * Rn * T^2;  % position from accel noise
L(4:6,1:3) = -Rn * T;          % velocity from accel noise
L(7:9,4:6) = -eye(3) * T;      % angle error from gyro noise

Qk = L * Sigma_z * L.';
Qk = (Qk + Qk.')/2;

end

function mhat = initLM(Y,mapInfo,Pcam_project,Ni)
K = size(Y,3);
mhat = zeros(3*K,1);
for k=1:K
    %find frames with kth LM visible
    Nk = sum(mapInfo.visibility(k,:));%number of frames with kth LM visible
    nk = find(mapInfo.visibility(k,:));%indices of frame for kth LM visible
    N = min(Ni,Nk);
    Yk = zeros(2,N);
    Pk = zeros(3,4,N);
    
    for id=1:N
        n = nk(id);
        Yk(:,id) = Y(:,n,k);  %select frames with kth LM
        Pk(:,:,id) = Pcam_project(:,:,n);
    end
    mhat(1+3*(k-1):3*k) = utils.TriangulateLS(Yk,Pk);
end
end

