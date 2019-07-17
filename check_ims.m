%% Start
addpath('blaser_util');
addpath('blaser_data/640_unified');
clear

n_val = 25;

%%
f = fopen('data.txt','r');
squareSize = 5; % mm
boardSize = [7,10];
threshold = 160;

figure;
for i = 1:n_val
    tr = fscanf(f, '[%f, %f, %f]\n', 3);
    rot = fscanf(f, '[%f, %f, %f, %f]\n', 4);
%     quat = [rot(4), rot(1), rot(2), rot(3)];
%     A{i} = [quat2rotm(quat), tr; 0,0,0,1];

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
    end
    
    laser_pixels = find_laser(I, threshold);
    if size(laser_pixels, 1) > 30
        disp("Laser line found!");
        plot(laser_pixels(:,1), laser_pixels(:,2), 'g.');
        drawnow;
    end
    
    pause(1.5);
end