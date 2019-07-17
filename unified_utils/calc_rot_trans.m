% function [rot, trans, J_R, J_t] = calc_rot_trans(state, H)
function [rot, trans] = calc_rot_trans(pts, world, state)
% Given the camera intrinsics and homography describing the point
%   projection, calculate the rotation and translation of the camera.

% K: camera intrinsic matrix (format [fx fy cx cy])
% H: Homography
% Returns:
%  rot- rotation of homography
%  trans- translation of homography
%  J_R: Jacobian of rotation. (9x12)
%       9 rows for rotation matrix
%       12 cols for 4 intrinsics + 8 homography
%  J_t: Jacobian of translation. (3x12)
%       3 rows for translation vector
%       12 cols for 4 intrinsics + 8 homography
%         omit h33 b/c its 1

H_vec = projective_transform(pts, world);
H = reshape([H_vec;1],3,3)';

fx = state(1);
fy = state(2);
cx = state(3);
cy = state(4);

Ki = [1/fx, 0,   -cx/fx;
       0,  1/fy, -cy/fy;
       0,   0,      1];
% This code copied from estimateCameraParameters.m
kh1 = Ki*H(:, 1);
kh2 = Ki*H(:, 2);
kh3 = Ki*H(:, 3);
lambda = 1 / norm(kh1);

% 3D rotation matrix
r1 = lambda * kh1;
r2 = lambda * kh2;
r3 = cross(r1, r2);
rot = [r1,r2,r3];
rot = rotationVectorToMatrix(vision.internal.calibration.rodriguesMatrixToVector(rot));

% translation vector
trans = lambda * kh3;


% % Jacobian of Ki * H1
% % Format (dims 3x12) 3 x [fx fy cx cy h11 h21 h31 h12 h22 h32 h13 h23]
% Jkh1 = zeros(3,12);
% Jkh1(1,1) = -kh1(1)/fx;
% Jkh1(2,2) = -kh1(2)/fy;
% Jkh1(1,3) = -kh1(3)/fx;
% Jkh1(2,4) = -kh1(3)/fy;
% Jkh1(:,5:7) = Ki;
% 
% % Jacobian of Ki * H2
% % Format (dims 3x7) 3 x [fx fy cx cy h11 h21 h31 h12 h22 h32 h13 h23]
% Jkh2 = zeros(3,12);
% Jkh2(1,1) = -kh2(1)/fx;
% Jkh2(2,2) = -kh2(2)/fy;
% Jkh2(1,3) = -kh2(3)/fx;
% Jkh2(2,4) = -kh2(3)/fy;
% Jkh2(:,8:10) = Ki;
% 
% % Jacobian of Ki * H3
% % Format (dims 3x6) 3 x [fx fy cx cy h11 h21 h31 h12 h22 h32 h13 h23]
% Jkh3 = zeros(3,12);
% Jkh3(1,1) = -kh3(1)/fx;
% Jkh3(2,2) = -kh3(2)/fy;
% Jkh3(1,3) = -kh3(3)/fx;
% Jkh3(2,4) = -kh3(3)/fy;
% Jkh3(1,11) = 1/fx;
% Jkh3(2,12) = 1/fy;
% 
% % Jacobian of lambda (gradient) 1x[fx fy cx cy h11 h21 h31 h12 h22 h32 h13 h23]
% Jl = -(lambda^3) * kh1'*Jkh1;
% 
% % Jacobian of r1-3+t (3x12)
% Jr1 = lambda * Jkh1 + kh1*Jl;
% Jr2 = lambda * Jkh2 + kh2*Jl;
% J_t =  lambda * Jkh3 + kh3*Jl;
% 
% % hodge1=[0,-kh1(3),kh1(2);kh1(3),0,-kh1(1);-kh1(2),kh1(1),0];
% % hodge2=[0,-kh2(3),kh2(2);kh2(3),0,-kh2(1);-kh2(2),kh2(1),0];
% hodge1=[0,-r1(3),r1(2);r1(3),0,-r1(1);-r1(2),r1(1),0];
% hodge2=[0,-r2(3),r2(2);r2(3),0,-r2(1);-r2(2),r2(1),0];
% Jr3 = lambda*(hodge1*Jkh2 - hodge2*Jkh1) + 2*r3*Jl/lambda;
% 
% % Jacobian of entire system now.
% J_R = [Jr1;Jr2;Jr3];
end
