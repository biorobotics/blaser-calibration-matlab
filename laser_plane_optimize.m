%% Start (all the things that change go here... i think.)
addpath('blaser_data/1280_verify')
addpath('blaser_util');
clear

load('1280_affine.mat')
load('a1001_1280_calib.mat')
laser_plane = [3.6048, -195.6298, -100, 3.6715319];
% laser_plane = [0, 0, -100, 0];

n_val = 12;

%% load file
f = fopen('data.txt','r');

A = cell(n_val, 1);
B = cell(n_val, 1);

imagePoints_all = [];
laserPoints_all = [];

for i = 1:n_val
    tr = fscanf(f, '[%f, %f, %f]\n', 3);
    rot = fscanf(f, '[%f, %f, %f, %f]\n', 4);
    quat = [rot(4), rot(1), rot(2), rot(3)];
    A{i} = [quat2rotm(quat), tr; 0,0,0,1];
    
    fname = fscanf(f, '%s\n', 1);
    B{i} = fname;
end

%% Interpret laser points
threshold = 150;
K = cameraParams.IntrinsicMatrix';

laser_points = cell(n_val, 1);

figure;
for i = 1:n_val
%     img = imread(B{i});
    img =  undistortImage(imread(B{i}), undistorter);
    [img, newOrigin] = undistortImage(img, cameraParams);
%     laser_pixels = extractPixelDataFromImg(img, threshold);
    laser_pixels = find_laser(img, 120);
    
    coeffs = polyfit(laser_pixels(:,1), laser_pixels(:,2), 1);

    dists = abs(polyval(coeffs, laser_pixels(:,1)) - laser_pixels(:,2));
    pts = laser_pixels(dists < 3, :);

    coeffs2 = polyfit(pts(:,1), pts(:,2), 1);
    dists = abs(polyval(coeffs2, laser_pixels(:,1)) - laser_pixels(:,2));
    pts = laser_pixels(dists < 1, :);
    
    clf;
%     imshow(img(:,:,1) > threshold);
    imshow(img);
    hold on;
    plot([0,1280], [0,1280]*coeffs2(1) + coeffs2(2), 'g')
    scatter(laser_pixels(:,1), laser_pixels(:,2), .5, 'b');
    scatter(pts(:,1), pts(:,2), 2, 'g');
    pause(1);

    laser_points{i} = pts;
%     pts3d = K \ [pts';ones(1,size(pts,1))];
%     z = -laser_plane(4) ./ (laser_plane(1:3) * pts3d);
%     laser_points{i} = pts3d .* z;
end


%% optimization?? idk setup step
threshold = 200;
K = cameraParams.IntrinsicMatrix';

scs = [];
fcs = [];
xyz = [];

for i = 1:n_val
    img = imread(B{i});
    [img, newOrigin] = undistortImage(img, cameraParams);
    laser_pixels = extractPixelDataFromImg(img, threshold);
    coeffs = polyfit(laser_pixels(:,1), laser_pixels(:,2), 1);
    
    dists = abs(polyval(coeffs, laser_pixels(:,1)) - laser_pixels(:,2));
    pts = laser_pixels(dists < 1, :);
    unitless = K \ [pts';ones(1,size(pts,1))];

    to_world = A{i}*aff;
    fc = (to_world(1:3,4) - checkerboard_loc')' * checkerboard_norm;
    sc = (to_world(1:3,1:3) * unitless)' * checkerboard_norm;
    
    scs = [scs;sc];
    fcs = [fcs;repmat(fc, size(sc))];
    xyz = [xyz; unitless'];

end

%% run iter
a = laser_plane(1);
b = laser_plane(2);
d = laser_plane(4);

prev_err = inf;
err = 1000000;
% err = norm(H);
% disp('Init error')
% disp(err)

mu = .8;
while prev_err/err-1 > 1E-9
    halp = xyz * [a;b;-100];
    H = -d./halp .* scs + fcs;
    J = [d*xyz(:,1) ./halp ./halp, d*xyz(:,2) ./halp ./halp, -1./halp] .*scs;
    
    jtj = J'*J;
    jth = J'*H;
    
    lhs = jtj + mu*eye(3);
    rhs = -jth;
    dX = double(lhs\rhs);
%     X = X+dX;
    a = a + dX(1);
    b = b + dX(2);
    d = d + dX(3);

    prev_err = err;
    err = norm(H);
    disp('Error:');
    disp(err);
end

disp('Done!');
disp('Final err:');
disp(err);
beep;

laser_plane = [a,b,-100,d]


%% Show everyting
stitched_cloud = [];

figure;
cm = hsv(n_val);
for i=1:n_val
    trplot(A{i}, 'length',.03, 'color', cm(i,:));
    hold on;
    
    camloc = A{i}*aff;
    plotCamera('Location', camloc(1:3,4), 'Orientation', camloc(1:3, 1:3)',...
               'Size', .01, 'Color', cm(i,:));
   
   pts3d = K \ [(laser_points{i})'; ones(1, size(laser_points{i},1))];
   z = -laser_plane(4) ./ (laser_plane(1:3) * pts3d);
   pts3d = pts3d .* z;
           
    warped = A{i} * aff * [pts3d; ones(1, size(pts3d, 2))];
    scatter3(warped(1,:), warped(2,:), warped(3,:), 2, cm(i,:));
    
    stitched_cloud = [stitched_cloud, warped(1:3,:)];
end
axis tight;
axis vis3d;

figure;
scatter3(stitched_cloud(1,:), stitched_cloud(2,:), stitched_cloud(3,:), 2);
axis equal;

%% evaluate goodness of fit
[fitresult, gof, plane] = fit3DFitPlane(stitched_cloud(1,:),stitched_cloud(2,:),stitched_cloud(3,:));
hold on;
axis tight;
axis(axis);

[x,y] = meshgrid(-1:0.02:1);
% x = 0.2*x; y=0.2*y;
z = (checkerboard_loc * checkerboard_norm - x*checkerboard_norm(1) - y*checkerboard_norm(2))/checkerboard_norm(3);
surf(x, y, z);

axis normal;


str = sprintf("[%.6f, %.6f, -100, %.6f]",...
    laser_plane(1), laser_plane(2), 1000*laser_plane(4));
disp("laser_plane:");
disp(str);
