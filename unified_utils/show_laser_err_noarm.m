function [] = show_laser_err_noarm(A,world,state)
%SHOW_LASER_ERR Display of the laser err in the 

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

laser_pts = [];
for i=1:n_im
    if size(A{i,3},1) == 0 || size(A{i,2},1) == 0
        continue
    end

    checkerboard = undistort_points(A{i,2}, state);
    [R, t] = calc_rot_trans(checkerboard, world, state);
    upts = undistort_points(A{i,3}, state);
    
    unitless = K \ [upts';ones(1,size(upts,1))]; %3xN vectors of points
    s = -D./(abc * unitless); %1xN vector of depths
    lpts = s.* unitless; %3xN vectors of points

    laser_pts = [laser_pts, R*lpts - R*t];
end

figure;
scatter3(laser_pts(1,:),laser_pts(2,:),laser_pts(3,:));
hold on;
axis tight;
axis(axis);

[x,y] = meshgrid(-1:0.02:1);
z = zeros(size(x));
surf(x, y, z);
axis normal;

end

