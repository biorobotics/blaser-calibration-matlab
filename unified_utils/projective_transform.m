% function [H, J_x, J_y] = projective_transform(pts, world)
function [H] = projective_transform(pts, world)
% pts are undistorted points (from image)
% world are world points (generated from checkerboard)
% Returns: H 8x1 (Homography transform)
%  J_x: Used to propogate derivative of Hij wrt each point
%  J_y: Used to propogate derivative of Hij wrt each point
% Each jacobian is 8xN, 8 components of 3x3 and N pts
% TODO: remove the for loop

npts = size(pts, 1);
x = world(:,1);
y = world(:,2);
v1 = ones(npts,1);
v0 = zeros(npts,1);
u = pts(:,1);
v = pts(:,2);

U = [u; v];
X = [x  y  v1 v0 v0 v0 -u.*x -u.*y;
     v0 v0 v0 x  y  v1 -v.*x -v.*y];

H = X\U;
% H = inv(X'*X)*X' * U;

% xtxi = pinv(X'*X);

% Mega chain rules lmao
% dxtx := d(xtx) = dX'*X + X'*dX
% dxtxi := d(inv(X'*X)) = -xtxi * dxtx * xtxi
% dH = (dxtxi*X' + xtxi*dX') * U + X\dU
% dH = ((-xtxi * (dX'*X + X'*dX) * xtxi)*X' + xtxi * dX') * U + X\dU;

% J_x = zeros(8, npts);
% J_y = zeros(8, npts);
% 
% for i=1:npts
%     dXx = zeros(2*npts, 8);
%     dXy = zeros(2*npts, 8);
%     
%     dXx (i,7:8) = -world(i,:);
%     dXy (npts+i,7:8) = -world(i,:);
%     
%     J_x(:,i) = (xtxi*dXx' - xtxi* (dXx'*X+X'*dXx) *xtxi*X') * U;
%     J_y(:,i) = (xtxi*dXy' - xtxi* (dXy'*X+X'*dXy) *xtxi*X') * U;
% end
% J_x = J_x + X\[eye(npts); zeros(npts)];
% J_y = J_y + X\[zeros(npts); eye(npts)];

end
