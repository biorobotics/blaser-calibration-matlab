function [] = show_laser_err(A,world,state)
%SHOW_LASER_ERR Display of the laser err in the 

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

[~, ch_loc, ch_ori] = handeye(A, world, state);
lpts = [];
for i=1:n_im
    if size(A{i,3},1) == 0
        continue
    end
    
    upts = undistort_points(A{i,3}, state);
    unitless = K \ [upts';ones(1,size(upts,1))]; %3xN vectors of points
    
    to_world = A{i,1}*aff;

    s = -D./(abc * unitless); %1xN vector of depths
    lpts = [lpts, s.* (to_world(1:3, 1:3)*unitless) + to_world(1:3, 4)]; %3xN vectors of points
end

figure;
scatter3(lpts(1,:),lpts(2,:),lpts(3,:));
hold on;
axis tight;
axis(axis);

[x,y] = meshgrid(-1:0.02:1);
z = (ch_loc' * ch_ori - x*ch_ori(1) - y*ch_ori(2))/ch_ori(3);
surf(x, y, z);
axis normal;

end

