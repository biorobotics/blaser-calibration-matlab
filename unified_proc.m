%% Start
clearvars -except blaserCalibDataPath
close all
addpath('unified_utils');

if exist('blaserCalibDataPath', 'var') ~= 1
    % Default calibration dataset
    global blaserCalibDataPath %#ok<TLEV>
    blaserCalibDataPath = './blaser_data/640_unified';
end

fprintf("Calibration dataset:\n    %s\n", blaserCalibDataPath);
n_val = numel(dir(fullfile(blaserCalibDataPath, '*.png'))) - 1;

fprintf("Dataset has %d images\n", n_val);

global blaserCalibDataPath %#ok<REDEFGG>

%% Load calibration data files
descFile = fullfile(blaserCalibDataPath, 'data.txt');
descFile = fopen(descFile,'r');
squareSize = .005; % meters
boardSize = [7,10];
[worldPoints] = generateCheckerboardPoints(boardSize,squareSize);
threshold = 160;

A = cell(n_val, 3);

figure;
for ii = 1:n_val
    tr = fscanf(descFile, '[%f, %f, %f]\n', 3);
    rot = fscanf(descFile, '[%f, %f, %f, %f]\n', 4);
    quat = [rot(4), rot(1), rot(2), rot(3)];
    A{ii, 1} = [quat2rotm(quat), tr; 0,0,0,1];

    fprintf('Image %d\n', ii);
    
    fname = fscanf(descFile, '%s\n', 1);
    I = imread(fullfile(blaserCalibDataPath, fname));
    
    clf;
    imshow(I);
    hold on;
    drawnow;
    
    [imagePoints,boardSize_detected] = detectCheckerboardPoints(I);
    if norm(boardSize_detected - boardSize) < 0.1
        disp("Checkerboard found!");
        plot(imagePoints(:,1), imagePoints(:,2), 'r');
        drawnow;
        
        A{ii,2} = imagePoints;
    end
    laser_pixels = find_laser_new(I);
    
    if numel(laser_pixels) > 1
        coeffs = polyfit(laser_pixels(:,1), laser_pixels(:,2), 1);
        dists = abs(polyval(coeffs, laser_pixels(:,1)) - laser_pixels(:,2));
        lpts = laser_pixels(dists < 2, :);

        coeffs2 = polyfit(lpts(:,1), lpts(:,2), 1);
        dists = abs(polyval(coeffs2, laser_pixels(:,1)) - laser_pixels(:,2));
        lpts = laser_pixels(dists < 1, :);

        if size(lpts, 1) > 30
            disp("Laser line found!");
            hold on
            % Plot found points
            plot(lpts(:,1), lpts(:,2), 'g.');
            % Plot interpolated line
            X = linspace(1, size(I, 2), size(I, 2));
            fline = polyval(coeffs, X);
            plot(X, fline, '--', 'color', [0 1 1])
            hold off
            drawnow;
            waitforbuttonpress;

            A{ii,3} = lpts;
        end
    end
end

%% Initialize camera intrinsics
imagePoints_all = [];
for ii=1:n_val
    if numel(A{ii,2}) ~= 0
    imagePoints_all = cat(3,imagePoints_all,A{ii,2});
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

while prev_err/err-1 > 1E-4
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

while prev_err/err-1 > 1E-5
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
fprintf('handeye calibration:\n<origin xyz="%.6f %.6f %.6f" rpy="%.6f %.6f %.6f"/>\n',...
    state(13),state(14),state(15),state(10),state(11),state(12));

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

while prev_err/err-1 > 1E-6
% while 1
    h = err_f(state);
    j = estimate_jac(err_f, state);

    ds = -pinv(j'*j + eye(size(j,2))) * j'*h;
    prev_err = err;
    state = state + ds;
    err = norm(err_f(state));
    fprintf('Til bound Error:%.9f\n',prev_err/err-1);  
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
    state(13),state(14),state(15),state(10),state(11),state(12));
