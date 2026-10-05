%camera parameters
fl0 = 3.5e-3;  %pixel focal length [m]
ps = 3e-6;
fl = fl0/ps;

asr = 4/3;
imageWidth = 6e-3;
imageHeight = imageWidth/asr;
imageWidthPx = imageWidth/ps;
imageHeightPx = imageHeight/ps;

L = [fl 0 0;0 fl 0;0 0 1];
Rbc = [1 0 0;0 0 1;0 -1 0];%body to camera (axis switch y to z)
phi = 2*atan(imageWidth/2/fl0);

Tc = 30e-3;
Fc = 1/Tc;

%imu parameters

sigma_a = 0.02;  % accel noise standard deviation per sample
sigma_g = 0.002; % gyro noise standard deviation per sample

Sigma_z = blkdiag(sigma_a^2*eye(3), ...
                  sigma_g^2*eye(3));
alpha = .001*ones(3,1);
gamma = [0.005;-0.003;0.01];

%sample interval [s]
Fs = 400; %Hz
T = 1/Fs;

