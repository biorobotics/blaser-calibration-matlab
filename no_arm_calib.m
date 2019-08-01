%% Start
addpath('unified_utils');
addpath('blaser_data/handheld_1280x720');
clear


%% Load calibration data files
n_val = count_lines('data.txt');
f = fopen('data.txt','r');
squareSize = .03; % meters
boardSize = [7,10];
[worldPoints] = generateCheckerboardPoints(boardSize,squareSize);
threshold = 130;

A = cell(n_val, 3);

figure;
for i = 1:n_val
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
    else
        continue
    end
       
    laser_pixels = find_laser(I, threshold); 
    if numel(laser_pixels) > 30
        disp("Laser line found!");
        plot(laser_pixels(:,1), laser_pixels(:,2), 'g.');
        drawnow;

        A{i,3} = laser_pixels;
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

while prev_err/err-1 > 1E-7
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

%% initial guess for laser points
[a,b,d] = init_laser_guess(A,worldPoints,state);
state = [state(1:9);a;b;1000*d];

%% Optimize laser plane and intrinsics all simultaneously
prev_err = inf;
err_f = @(st) laser_plane_err_noarm(A, worldPoints, st);

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
    fprintf('Error: %.7f\n', err);
end
disp("DONE");
fprintf('Final err: %.7f %.7f %.7f\n', laser_plane_err_noarm(A, worldPoints, state,true));

%% Show laser error
show_laser_err_noarm(A,worldPoints,state);
title('Laser plane error (intrinsic + laser)');

show_reproj_err(A,worldPoints,state);
title('Reprojection error (intrinsic + laser)');

%% Display results in a nice way
fprintf('Intransic matrix:\n[%.5f, 0, %.5f, 0, %.5f, %.5f, 0, 0, 1]\n',...
    state(1), state(3), state(2), state(4)); % fx, cx, fy, cy
fprintf('Distortion coefficients:\n[%.5f, %.5f, %.5f, %.5f, %.5f]\n',...
    state(5),state(6),state(7),state(8),state(9)); % k1 k2 p1 p2 k3
fprintf('Laser Plane:\n[%.5f, %.5f, -100, %.5f]\n',...
    state(10), state(11), state(12)); % a b d


