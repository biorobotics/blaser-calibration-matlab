%% Start
addpath('unified_utils');
addpath('blaser_data/1280_unified');
clear

n_val = 30;

%% Load calibration data files
f = fopen('data.txt','r');
squareSize = .005; % meters
boardSize = [7,10];
[worldPoints] = generateCheckerboardPoints(boardSize,squareSize);
threshold = 160;

A = cell(n_val, 3);

figure;
for i = 1:n_val
    tr = fscanf(f, '[%f, %f, %f]\n', 3);
    rot = fscanf(f, '[%f, %f, %f, %f]\n', 4);
    quat = [rot(4), rot(1), rot(2), rot(3)];
    A{i, 1} = [quat2rotm(quat), tr; 0,0,0,1];

    fprintf('Image %d\n', i);
    
    fname = fscanf(f, '%s\n', 1);
    I = imread(fname);
    
    clf;
    imshow(I);
    hold on;
    drawnow;
    
    [imagePoints,boardSize_detected] = detectCheckerboardPoints(I);
    if norm(boardSize_detected - boardSize) < 0.1
        disp("Checkerboard found!");
        plot(imagePoints(:,1), imagePoints(:,2), 'r');
        drawnow;
        
        A{i,2} = imagePoints;
    end
    
    laser_pixels = find_laser(I, threshold);
    
    if numel(laser_pixels) > 1
        coeffs = polyfit(laser_pixels(:,1), laser_pixels(:,2), 1);
        dists = abs(polyval(coeffs, laser_pixels(:,1)) - laser_pixels(:,2));
        lpts = laser_pixels(dists < 3, :);

        coeffs2 = polyfit(lpts(:,1), lpts(:,2), 1);
        dists = abs(polyval(coeffs2, laser_pixels(:,1)) - laser_pixels(:,2));
        lpts = laser_pixels(dists < 1, :);

        if size(lpts, 1) > 30
            disp("Laser line found!");
            plot(lpts(:,1), lpts(:,2), 'g.');
            drawnow;

            A{i,3} = lpts;
        end
    end
end

%% Initialize camera intrinsics
imagePoints_all = [];
for i=1:n_val
    if numel(A{i,2}) ~= 0
    imagePoints_all = cat(3,imagePoints_all,A{i,2});
    end
end
cp = estimateCameraParameters(imagePoints_all, worldPoints,...
    'NumRadialDistortionCoefficients',3,'EstimateTangentialDistortion',true);

state = [cp.FocalLength  cp.PrincipalPoint cp.RadialDistortion(1:2),...
    cp.TangentialDistortion cp.RadialDistortion(3)]';

%% Optimize camera intrinsics by minimizing reprojection error
prev_err = inf;
err_f = @(st) vectorize_reproj(A, worldPoints, st);

err = norm(err_f(state));
fprintf('Initial error: %.7f\n', err);

while prev_err/err-1 > 1E-5
    h = err_f(state);
    j = estimate_jac(err_f, state);

    ds = -pinv(j'*j + eye(size(j'*j))) * j'*h;
    prev_err = err;
    state = state + ds;
    err = norm(err_f(state));
    fprintf('Error: %.7f\n', err);
end
disp("DONE");
fprintf('Final err: %.7f\n', vectorize_reproj(A,worldPoints,state,true));

%% Show reprojection errors of checkerboard
show_reproj_err(A,worldPoints,state);
title('Reprojection error (only camera)');

%% init hand-eye
state = [state(1:9)', 0,0,0, 0,0,0]';

%% Optimize handeye simultaneously with intrinsics
prev_err = inf;
err_f = @(st) handeye(A, worldPoints, st);

err = norm(err_f(state));
fprintf('Initial error: %.7f\n', err);

while prev_err/err-1 > 1E-4
    h = err_f(state);
    j = estimate_jac(err_f, state);
   
    ds = -pinv(j'*j + eye(size(j,2))) * j'*h;
    prev_err = err;
    state = state + ds;
    err = norm(err_f(state));
    fprintf('Error: %.7f\n', err);
end
disp("DONE");
fprintf('Final err: %.7f %.7f\n', handeye(A, worldPoints, state, true));

%% Show hand-eye fit
show_handeye(A, worldPoints, state);
title('Handeye (instrinsic + handeye)');

show_reproj_err(A,worldPoints,state);
title('Reprojection error (intrinsic + handeye)');

%% Initialize laser plane params
state = [state(1:15); 4.099360;-196.800802;3613.930297];

%% Optimize laser plane and handeye and intrinsics all simultaneously
prev_err = inf;
err_f = @(st) laser_plane_err(A, worldPoints, st);

err = norm(err_f(state));
fprintf('Initial error: %.7f\n', err);

while prev_err/err-1 > 1E-8
% while 1
    h = err_f(state);
    j = estimate_jac(err_f, state);

    ds = -pinv(j'*j + eye(size(j,2))) * j'*h;
    prev_err = err;
    state = state + ds;
    err = norm(err_f(state));
    fprintf('Error: %.7f %.7f %.7f %.7f\n', laser_plane_err(A, worldPoints, state,true));
end
disp("DONE");
fprintf('Final err: %.7f %.7f %.7f %.7f\n', laser_plane_err(A, worldPoints, state,true));

%% Show laser error
show_laser_err(A, worldPoints, state);
title('Laser plane error (instrinsic + handeye + laser)');

show_handeye(A, worldPoints, state);
title('Handeye (instrinsic + handeye + laser)');

show_reproj_err(A,worldPoints,state);
title('Reprojection error (intrinsic + handeye + laser)');

%% Display results in a nice way
fprintf('Intransic matrix:\n[%.5f, 0, %.5f, 0, %.5f, %.5f, 0, 0, 1]\n',...
    state(1), state(3), state(2), state(4)); % fx, cx, fy, cy
fprintf('Distortion coefficients:\n[%.5f, %.5f, %.5f, %.5f, %.5f]\n',...
    state(5),state(6),state(7),state(8),state(9)); % k1 k2 p1 p2 k3
fprintf('Laser Plane:\n[%.5f, %.5f, -100, %.5f]\n',...
    state(16), state(17), state(18)); % a b d

fprintf('handeye calibration:\n<origin xyz="%.6f %.6f %.6f" rpy="%.6f %.6f %.6f"/>\n',...
    state(13),state(14),state(15),state(12),state(11),state(10));
