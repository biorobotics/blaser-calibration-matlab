function [] = unified_proc(calibDataPath, squareSize, boardSize)
%UNIFIED_PROC Unified calibration process for camera parameters,
%             laser plane, handeye transform based on optimization.
%             Note that running with function call is discoraged,
%             best practice is to run section by section.
%% Start
close all
addpath('unified_utils');

if ~exist('calibDataPath', 'var')
    % Default calibration dataset
    calibDataPath = './blaser_data/640_unified';
end

if ~exist('squareSize', 'var')
    squareSize = .005; % meters
end

if ~exist('boardSize', 'var')
    boardSize = [7, 10];
end

minLaserStripePts = 70;

f_listing = dir(fullfile(calibDataPath, '*.png'));
n_val = numel(f_listing);

if n_val < 5
    fprintf("Warning: dataset is too small, exiting...\n");
    return
end

fprintf("------------------------\n");
fprintf("Calibration dataset:\n    %s\n", calibDataPath);
fprintf("Dataset size:\n    %d\n", n_val);
fprintf("Checkerboard size:\n    %.4f\n", squareSize);
fprintf("Checkerboard dims:\n    %dx%d\n", boardSize(1), boardSize(2));

%% Load calibration data files
descFile = fullfile(calibDataPath, 'data.txt');
descFile = fopen(descFile,'r');

handEyeEnabled = false;
if descFile > 0
    handEyeEnabled = true;
end
fprintf("Handeye enabled:\n    %d\n", handEyeEnabled);

[worldPoints] = generateCheckerboardPoints(boardSize,squareSize);

A = cell(n_val, 3);

imgnames = repmat({''}, 1, n_val);

for ii = 1:n_val
    % Get image files and hand pose if applicable
    if handEyeEnabled
        tr = fscanf(descFile, '[%f, %f, %f]\n', 3);
        rot = fscanf(descFile, '[%f, %f, %f, %f]\n', 4);
        quat = [rot(4), rot(1), rot(2), rot(3)];
        A{ii, 1} = [quat2rotm(quat), tr; 0,0,0,1];
        imgnames{ii} = fscanf(descFile, '%s\n', 1);
    else
        imgnames{ii} = f_listing(ii).name;
    end
    
    % Get checkerboard points
    fname = imgnames{ii};
    J = imread(fullfile(calibDataPath, fname));
    
    [imagePoints,boardSize_detected] = detectCheckerboardPoints(J);
    if norm(boardSize_detected - boardSize) < 0.1
        disp("Checkerboard found!");
        
        
        A{ii,2} = imagePoints;
    end
end

%% Initialize camera intrinsics
imagePoints_all = [];
for ii=1:n_val
    if numel(A{ii,2}) ~= 0
    imagePoints_all = cat(3,imagePoints_all,A{ii,2});
    end
end
cp = estimateCameraParameters(imagePoints_all, worldPoints, ...
    'NumRadialDistortionCoefficients', 3, ...
    'EstimateTangentialDistortion', true);

state = [cp.FocalLength, cp.PrincipalPoint, cp.RadialDistortion(1:2), ...
         cp.TangentialDistortion, cp.RadialDistortion(3)]';

%% Get laser stripes
figure;
for ii = 1:n_val
    fname = imgnames{ii};
    
    % Undistort image before processing
    I = imread(fullfile(calibDataPath, fname));
    [J, newOrigin] = undistortImage(I, cp, 'outputView', 'full');
    
    clf;
    imshow(J);
    hold on;
    drawnow;
    
    imagePoints = A{ii,2};
    if numel(imagePoints) ~= 0
        imagePoints = undistortPoints(imagePoints, cp);
        plot(imagePoints(:,1) - newOrigin(1), ...
            imagePoints(:,2) - newOrigin(2), 'g*-');
        drawnow;
    end
    
    % Find laser lines
    laser_pixels = find_laser_new(J);
    
    if numel(laser_pixels) > 1
        coeffs = polyfit(laser_pixels(:,1), laser_pixels(:,2), 1);
        dists = abs(polyval(coeffs, laser_pixels(:,1)) - laser_pixels(:,2));
        lpts = laser_pixels(dists < 2, :);

        coeffs2 = polyfit(lpts(:,1), lpts(:,2), 1);
        dists = abs(polyval(coeffs2, laser_pixels(:,1)) - laser_pixels(:,2));
        lpts = laser_pixels(dists < 1, :);

        if size(lpts, 1) > minLaserStripePts
            disp("Laser line found!");
            hold on
            % Plot found points
            plot(lpts(:,1), lpts(:,2), 'g.');
            % Plot interpolated line
            X = linspace(1, size(J, 2), size(J, 2));
            fline = polyval(coeffs, X);
            plot(X, fline, '--', 'color', [0 1 1])
            hold off
            drawnow;
            waitforbuttonpress;

            A{ii,3} = lpts;
        end
    end
end

if descFile > 0
    fclose(descFile);
end



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
if handEyeEnabled
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
    fprintf('Final err: %.7f %.7f\n', ...
        handeye(A, worldPoints, state, true));
end

%% Show hand-eye fit
if handEyeEnabled
    fprintf('handeye calibration:\n<origin xyz="%.6f %.6f %.6f" rpy="%.6f %.6f %.6f"/>\n',...
        state(13),state(14),state(15),state(10),state(11),state(12));
    
    show_handeye(A, worldPoints, state);
    title('Handeye (instrinsic + handeye)');
    
    show_reproj_err(A, worldPoints, state);
    title('Reprojection error (intrinsic + handeye)');
end

%% Initialize laser plane params
state = [state(1:15); 4.099360;-196.800802;3613.930297];

% [a,b,d] = init_laser_guess(A,worldPoints,state);
% state = [state(1:15);a;b;1000*d];

%% Optimize laser plane and handeye and intrinsics all simultaneously
prev_err = inf;
if handEyeEnabled
    err_f = @(st) laser_plane_err(A, worldPoints, st);
else
    err_f = @(st) laser_plane_err_noarm(A, worldPoints, st);
end

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
    if handEyeEnabled
        fprintf('Error: %.7f %.7f %.7f %.7f\n', ...
            laser_plane_err(A, worldPoints, state,true));
    else
        fprintf('Error: %.7f %.7f %.7f\n', ...
            laser_plane_err_noarm(A, worldPoints, state,true));
    end
end
disp("DONE");
if handEyeEnabled
    fprintf('Final err: %.7f %.7f %.7f %.7f\n', ...
        laser_plane_err(A, worldPoints, state,true));
else
    fprintf('Final err: %.7f %.7f %.7f\n', ...
        laser_plane_err_noarm(A, worldPoints, state,true));
end

%% Show laser error
if handEyeEnabled
    show_laser_err(A, worldPoints, state);
else
    show_laser_err_noarm(A, worldPoints, state);
end
title('Laser plane error (instrinsic + handeye + laser)');

if handEyeEnabled
    show_handeye(A, worldPoints, state);
    title('Handeye (instrinsic + handeye + laser)');
end

show_reproj_err(A,worldPoints,state);
title('Reprojection error (intrinsic + handeye + laser)');

%% Display results in a nice way
fprintf('Intrinsic matrix:\n[%.5f, 0, %.5f, 0, %.5f, %.5f, 0, 0, 1]\n',...
    state(1), state(3), state(2), state(4)); % fx, cx, fy, cy
fprintf('Distortion coefficients:\n[%.5f, %.5f, %.5f, %.5f, %.5f]\n',...
    state(5),state(6),state(7),state(8),state(9)); % k1 k2 p1 p2 k3
fprintf('Laser Plane:\n[%.5f, %.5f, -100, %.5f]\n',...
    state(16), state(17), state(18)); % a b d

if handEyeEnabled
    fprintf('handeye calibration:\n<origin xyz="%.6f %.6f %.6f" rpy="%.6f %.6f %.6f"/>\n',...
        state(13),state(14),state(15),state(10),state(11),state(12));
end

end