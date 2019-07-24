
function [errs, R, t] = reproj_err(pts, world, state)
% From the pts, world, and state, calculate the current reprojection errors
% Will undistort the points first.
% Project world points into camera frame using intrinsic matrix
%   & Homography (from pts to world)

% pts- Nx2 vector of points
% world- Nx2 'ground truth' values
% state- 9x1 vector of parameters [fx fy cx cy k1 k2 p1 p2 k3]
% Returns:
%  err- 2Nx1 vector of reprojection errors (x and y)
%  jac- 2Nx9 matrix (derivative of err by state)

fx = state(1);
fy = state(2);
cx = state(3);
cy = state(4);
% k1 = state(5);
% k2 = state(6);
% p1 = state(7);
% p2 = state(8);
% k3 = state(9);

npts = size(pts, 1);

% [upts, J_uptx, J_upty] = undistort_points(pts, [640, 480], cx, cy, k1, k2, p1, p2, k3);
% [H_vec, J_Hx, J_Hy] = projective_transform(upts, world);
upts = undistort_points(pts, state);

[R, t] = calc_rot_trans(upts, world, state);

K = [fx,  0, 0;
     0,  fy, 0;
     cx, cy, 1];

proj_3 = [world, zeros(npts,1), ones(npts,1)] * [R; t'] * K; % Nx3
proj = proj_3(:,1:2) ./ proj_3(:,3);

diffs = proj - upts;
errs = [diffs(:,1); diffs(:,2)]; % 2Nx1

% jac = zeros(2*npts, 9);
% % undistort contribution
% jac(1:npts, 3:9) = -J_uptx;
% jac(npts+1:2*npts, 3:9) = -J_upty;
% 
% % proj contribution?
% J_proj3 = [J_R(7:9,:); J_t(3,:)];
% J_proj1 = [R;t'] * [[1,0,0,0;0,0,0,0;0,0,1,0], zeros(3,8)]; % 3x12 Jacobian of K(:,1)
end
