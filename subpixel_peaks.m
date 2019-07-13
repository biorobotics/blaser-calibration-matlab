%% Start
addpath('blaser_util');
addpath('blaser_data/1280_verify');
clear

threshold = 150;

%%
% v = VideoReader('a1001_1280x960_calib_final.mp4');
% n_frame = round(v.Duration*v.FrameRate)-2;

% i = 1;
I = imread('im0.png');
[pts] = find_laser(I, 150);
% disp('wow slow');
figure;
imshow(I);

laser_pixels = extractPixelDataFromImg(I, threshold);
coeffs = polyfit(laser_pixels(:,1), laser_pixels(:,2), 1);

dists = abs(polyval(coeffs, laser_pixels(:,1)) - laser_pixels(:,2));
lpts = laser_pixels(dists < 1, :);

hold on
scatter(pts(:,1), pts(:,2), 1, 'g.');

%%
figure;
imshow(I);
hold on;
scatter(lpts(:,1), lpts(:,2), 1, 'g.')
% scatter(pts(:,1), pts(:,2), 1, 'g.')
