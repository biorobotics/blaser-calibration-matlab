function [resid] = laser_plane_err(A, world, state)
%LASER_PLANE Gets laser plane error from state laser plane parameters
%   Takes checkerboard location and orientation from handeye.
%   Error is deviation from plane; also fit a line to 3d points and
%    get laser line direction & take error in orientation
% State- 18x1 vector of values to optimize over
%  [fx fy cx cy k1 k2 p1 p2 k3 r p w x y z A B D] (C=-100)

n_im = size(A,1);
fx = state(1);
fy = state(2);
cx = state(3);
cy = state(4);
r = state(10);
p = state(11);
w = state(12);
x = state(13);
y = state(14);
z = state(15);
abc = [state(16), state(17), -100];
D = state(18)/1000;
aff = [eul2rotm([w,p,r]),[x;y;z];0,0,0,1];
K = [fx,  0, cx;
      0, fy, cy;
      0,  0, 1];

  
[resid, ch_loc, ch_ori] = handeye(A, world, state);

for i=1:n_im
    if size(A{i,3},1) == 0
        continue
    end
    
    upts = undistort_points(A{i,3}, state);
%     upts = A{i,3};
    unitless = K \ [upts';ones(1,size(upts,1))]; %3xN vectors of points
    
    to_world = A{i,1}*aff;

    s = -D./(abc * unitless); %1xN vector of depths
    lpts = s.* (to_world(1:3, 1:3)*unitless) + to_world(1:3, 4); %3xN vectors of points
    
    % out of plane error
    errs = (lpts - ch_loc)' * ch_ori; % Nx1 residual vector
    
    % line direction error
    l_mean = mean(lpts, 2);
    dX = lpts - l_mean;
    C = (dX*dX')/(size(lpts,1)-1);
    [R,~]=svd(C,0);
    dir = R(:,1);
    % dir' * ch_ori;

    resid = [resid;20*errs; 50* dir' * ch_ori];
end

end

