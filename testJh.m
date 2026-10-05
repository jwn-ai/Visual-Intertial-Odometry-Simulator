function [H_numeric,H_single] = testJh(x,y0,rk,fl,Rbc)
        epsilon = 1e-6; % Tiny perturbation angle
H_numeric = zeros(2,15);
% 1. Perturb Position (Columns 1:3) and Velocity (Columns 4:6)
    Rn = utils.q2R(x(7:10));

for i = 1:6
    x_pert = x;
    x_pert(i) = x_pert(i) + epsilon;
    y = h(rk,Rn*Rbc,x_pert(1:3),fl);

    H_numeric(:, i) = (y - y0) / epsilon;
end

% 2. Perturb 3D Orientation Error (Columns 7:9)
for i = 1:3
    dtheta = zeros(3, 1);
    dtheta(i) = epsilon;
    
    % Map the 3D angle change to a local 4D delta-quaternion
    dq = [1; 0.5 * dtheta(1); 0.5 * dtheta(2); 0.5 * dtheta(3)];
    dq = dq / norm(dq);
    
    x_pert = x;
x_pert(7:10) = correctq(dtheta, x(7:10));
Rn_pert = utils.q2R(x_pert(7:10));
y = h(rk, Rn_pert*Rbc, x_pert(1:3), fl);

    H_numeric(:, 6 + i) = (y - y0) / epsilon;
end

% 3. Perturb Gyro Bias (Columns 10:12) and Accel Bias (Columns 13:15)
% Note the offset: state index is 11:16, but Jacobian index is 10:15
for i = 1:6
    state_idx = 10 + i;     % Maps to x(11) through x(16)
    jac_idx = 9 + i;        % Maps to H(:, 10) through H(:, 15)
    
    x_pert = x;
    x_pert(state_idx) = x_pert(state_idx) + epsilon;
    y = h(rk,Rn*Rbc,x_pert(1:3),fl);

    H_numeric(:, jac_idx) = (y - y0) / epsilon;
end



% Slot 1: Perturb Vehicle Pose (Columns 1:15) -> Use our working code blocks
% ... [Your working 1:15 vehicle perturbation loop here] ...

% Slot 2: Perturb the EXACT Landmark State Coordinates (Columns for Landmark i)

H_single = zeros(2,3);

for j = 1:3
    rk_pert = rk;
    rk_pert(j) = rk_pert(j) + epsilon;

    z_pert = h(rk_pert, Rn*Rbc, x(1:3), fl);
    H_single(:,j) = (z_pert - y0) / epsilon;
end
end

function hx = h(rk,R,r,fl)
    pk = utils.TF(rk,R,r);
    hx = utils.PIz(pk,fl);
end
