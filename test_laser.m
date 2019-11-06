function [] = test_laser(calibDataPath)

%% Load dataset
close all
addpath('unified_utils');

if ~exist('calibDataPath', 'var')
    % Default calibration dataset
    calibDataPath = './blaser_data/640_unified';
end

f_listing = dir(fullfile(calibDataPath, '*.png'));
n_val = numel(f_listing);
imgnames = repmat({''}, 1, n_val);
for ii = 1:n_val
    imgnames{ii} = f_listing(ii).name;
end

minLaserStripePts = 70;

A = cell(n_val, 3);

fprintf("------------------------\n");
fprintf("Calibration dataset:\n    %s\n", calibDataPath);
fprintf("Dataset size:\n    %d\n", n_val);

% Note: this test assumes image is undistorted!
% This test is intended for testing stripe extraction pipeline only.

%% Laser test
figure;
for ii = 1:n_val
      fname = imgnames{ii};
      fprintf('Image %s\n', fname);
      I = imread(fullfile(calibDataPath, fname));

      imshow(I);
      hold on;
      drawnow;

    laser_pixels = find_laser_new(I);
    plot(laser_pixels(:,1), laser_pixels(:,2), 'r.');

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
            X = linspace(1, size(I, 2), size(I, 2));
            fline = polyval(coeffs, X);
            plot(X, fline, '--', 'color', [0 1 1])
            hold off
            drawnow;

            A{ii,3} = lpts;
        end

        waitforbuttonpress;
    end
        
end
