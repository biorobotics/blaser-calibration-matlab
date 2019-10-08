%% Start
addpath('unified_utils');
addpath('blaser_data/640_unified');
clear

n_val =2;
%% Load calibration data files
f = fopen('data.txt','r');
squareSize = .004; % meters
boardSize = [7,10];
% [worldPoints] = generateCheckerboardPoints(boardSize,squareSize);
threshold = 160;

A = cell(n_val, 3);

figure;
for i = n_val:n_val
%     tr = fscanf(f, '[%f, %f, %f]\n', 3);
%     rot = fscanf(f, '[%f, %f, %f, %f]\n', 4);
%     quat = [rot(4), rot(1), rot(2), rot(3)];
%     A{i, 1} = [quat2rotm(quat), tr; 0,0,0,1];
% 
      fprintf('Image %d\n', i);
%     
%       fname = fscanf(f, '%s\n', 1)
      fname = append("im", int2str(i),".png");
      I = imread(fname);
%     
%     clf;
      imshow(I);
      hold on;
      drawnow;
%     
%     [imagePoints,boardSize_detected] = detectCheckerboardPoints(I);
%     if norm(boardSize_detected - boardSize) < 0.1
%         disp("Checkerboard found!");
%         plot(imagePoints(:,1), imagePoints(:,2), 'r');
%         drawnow;
%         
%         A{i,2} = imagePoints;
%     end
    
    laser_pixels = find_laser_new(I);
    disp(size(laser_pixels));
    
    if numel(laser_pixels) > 1
%         coeffs = polyfit(laser_pixels(:,1), laser_pixels(:,2), 1);
%         dists = abs(polyval(coeffs, laser_pixels(:,1)) - laser_pixels(:,2));
%         lpts = laser_pixels(dists < 100, :);
% 
%         coeffs2 = polyfit(lpts(:,1), lpts(:,2), 1);
%         dists = abs(polyval(coeffs2, laser_pixels(:,1)) - laser_pixels(:,2));
        lpts = laser_pixels;
        disp(size(lpts));
        if size(lpts, 1) > 30
            disp("Laser line found!");
            plot(lpts(:,1), lpts(:,2), 'g.');
            drawnow;

            A{i,3} = lpts;
        end
    end
end
