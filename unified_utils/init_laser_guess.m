function [a,b,d] = init_laser_guess(A, world, state)
%INIT_LASER_GUESS Summary of this function goes here
%   Detailed explanation goes here

n_im = size(A,1);
fx = state(1);
fy = state(2);
cx = state(3);
cy = state(4);

% abc = [state(16), state(17), -100];
% D = state(18)/1000;
% aff = [eul2rotm([w,p,r]),[x;y;z];0,0,0,1];
K = [fx,  0, cx;
      0, fy, cy;
      0,  0, 1];

sand = [];
for i=1:n_im
    if size(A{i,3},1) == 0
        continue
    end
    
    checkerboard = undistort_points(A{i,2}, state);
    %lpts = undistort_points(A{i,3}, state);
    lpts = A{i,3};
    
    [R, t] = calc_rot_trans(checkerboard, world, state);
    R=R';
    d = -t'*R(:,3);
    
    point_3d = K\[lpts';ones(1,size(lpts,1))];
    z = -d./(point_3d'*R(:,3));
        
    sand = [sand, point_3d.*[z';z';z']];
end

mpt = mean(sand,2);
[eigv,~,~] = svd(cov(sand'));

norm = eigv(:,3);

% Scale 'c' to -100
sf = -100/norm(3);
a = norm(1)*sf;
b = norm(2)*sf;
d = -norm'*mpt *sf;

end

