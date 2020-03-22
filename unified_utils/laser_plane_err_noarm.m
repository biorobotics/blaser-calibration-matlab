function [resid] = laser_plane_err_noarm(A, world, state, summarize)
%LASER_PLANE Gets laser plane error from state laser plane parameters
%   Assumes checkerboard location and orientation as identity.
%   Error is deviation from plane; also fit a line to 3d points and
%    get laser line direction & take error in orientation
% State- 18x1 vector of values to optimize over
%  [fx fy cx cy k1 k2 p1 p2 k3 r p w x y z A B D] (C=-100)

if ~exist('summarize','var')
    summarize = false;
end

n_im = size(A,1);
fx = state(1);
fy = state(2);
cx = state(3);
cy = state(4);
abc = [state(16), state(17), -100];
D = state(18)/1000;
K = [fx,  0, cx;
      0, fy, cy;
      0,  0, 1];

rpr_resid = [];
plane_resid = [];
dir_resid = [];
for i=1:n_im
    if size(A{i,3},1) == 0 || size(A{i,2},1) == 0
        continue
    end

    [rpr_err, R, t] = reproj_err(A{i,2}, world, state);
    rpr_resid = [rpr_resid;.1*rpr_err];
    %upts = undistort_points(A{i,3}, state);
    upts = A{i,3};
    
    unitless = K \ [upts';ones(1,size(upts,1))]; %3xN vectors of points
    s = -D./(abc * unitless); %1xN vector of depths
    lpts = s.* unitless; %3xN vectors of points

    lpts = R*lpts - R*t;
%     B = [R,-R*t; 0,0,0,1];
%     to_world = inv(B);

    % line direction error
    [R,~]=svd(cov(lpts'),0);
    dir = R(:,1);


%     lpz = [lpz,R*lpts - R*t];
    plane_resid = [plane_resid; 50*lpts(3,:)']; % Nx1 residual vector
    dir_resid = [dir_resid; 20*dir(3)];
end

if summarize
    resid = [norm(rpr_resid);norm(plane_resid);norm(dir_resid)];
else
    resid = [rpr_resid;plane_resid;dir_resid];
end

end

