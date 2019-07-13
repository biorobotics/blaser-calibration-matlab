%% Start
addpath('blaser_util');
addpath('blaser_data');
clear

%% Load images for camera calibration (just to get undistort coeffs)
v = VideoReader('a1001_1280x960_calib_final.mp4');
% v = VideoReader('a1001_640x480_calib_final.mp4');

squareSize = 5; % mm
skip_factor = 10; % integer plz
boardSize = [7,10];
[worldPoints] = generateCheckerboardPoints(boardSize,squareSize);
n_frame = round(v.Duration*v.FrameRate)-2;
disp(n_frame);
threshold = 200;
imagePoints_all = [];
frame = [];
i = 1;

while i < n_frame
    fprintf('image %d\n', i)
    I = read(v,i);
%     I = read(v,i);
    clf;
    imshow(I);
    drawnow;

    % get checkerboard
    [imagePoints,boardSize_detected] = detectCheckerboardPoints(I);

    if norm(boardSize - boardSize_detected)>0.1
        i = i+1;
        continue
    end
    
    hold on;
    plot(imagePoints(:,1), imagePoints(:,2), 'r');
    drawnow;

    scatter(imagePoints(1,1), imagePoints(1,2), 'p');
    drawnow;
    
    imagePoints_all = cat(3,imagePoints_all,flip(imagePoints,1)); %why flip
    frame = [frame,i];
    i = i+skip_factor;
end

%%
[undistorter,~,~] = estimateCameraParameters(imagePoints_all,worldPoints, ...
                                                'NumRadialDistortionCoefficients', 3,  'EstimateTangentialDistortion', true);


%% Load all the images again, but undistort first

imagePoints_undistort = [];
del_frame = [];
figure;
for i = 1:size(frame, 2)
    disp(frame(i))
    I = read(v,frame(i));
    J = undistortImage(I, undistorter);
    
    clf;
    imshow(J);
    hold on;
    
    % get checkerboard
    [imagePoints,boardSize_detected] = detectCheckerboardPoints(J);
    if norm(boardSize - boardSize_detected)>0.1
        del_frame = [del_frame, i];
        continue
    end
    
    hold on;
    plot(imagePoints(:,1), imagePoints(:,2), 'r');
    drawnow;
    
    imagePoints_undistort = cat(3,imagePoints_undistort,flip(imagePoints,1)); %why flip
end
frame(del_frame) = [];

%%
[cameraParams,imagesUsed,estimationErrors] = estimateCameraParameters(imagePoints_undistort,worldPoints, ...
                                                'NumRadialDistortionCoefficients', 3,  'EstimateTangentialDistortion', true);
%
%% something i gues
figure;
subplot(3,2,[1 3]); showExtrinsics(cameraParams, 'CameraCentric');
subplot(3,2,[2 4]); showExtrinsics(cameraParams, 'PatternCentric');
subplot(3,2,[5 6]);showReprojectionErrors(cameraParams);
% displayErrors(estimationErrors, cameraParams);


%% get reprojection errors
errs = [];
c_by_im = [];
c_by_pt = [];

cm = hsv(size(frame,2));
% cm2 = hsv(54);
for i = 1:size(frame, 2)
    disp(frame(i))
    I = read(v,frame(i));
    J = undistortImage(I, undistorter);
    
%     clf;
%     imshow(J);
%     hold on;
%     plot(imagePoints_undistort(:,1,i), imagePoints_undistort(:,2,i), 'go');
    warped_pts = [worldPoints, zeros(size(worldPoints,1),1), ones(size(worldPoints,1),1)]...
        * [cameraParams.RotationMatrices(:,:,i);cameraParams.TranslationVectors(i,:)]...
        * cameraParams.IntrinsicMatrix;
    warped_pts = warped_pts ./ warped_pts(:, 3);
    errs = [errs;warped_pts(:,1:2) - imagePoints_undistort(:,:,i)];
    c_by_im = [c_by_im;repmat(cm(i,:), size(warped_pts,1), 1)];
    c_by_pt = [c_by_pt;hsv(size(warped_pts,1))];
%     scatter(warped_pts(:,1), warped_pts(:,2), 'r+');
%     pause(0.5);
end


%%
figure;
subplot(1,2,1);
scatter(errs(:,1), errs(:,2), 2, c_by_im);
hold on;
merrs = [];
for i=1:size(frame,2)
    merrs = [merrs; mean(errs(54*i-53:54*i,:),1)];
end
scatter(merrs(:,1), merrs(:,2), 25, hsv(size(frame,2)), '*');

subplot(1,2,2);
scatter(errs(:,1), errs(:,2), 2, c_by_pt);
hold on;
merrs = mean(reshape(errs, 54, size(frame,2), 2), 2);
scatter(merrs(:,:,1), merrs(:,:,2), 25, hsv(54), '*');


%% print params

disp(" ");
disp("Intrinsic matrix:");
str = "[";
for i = 1:8
    if undistorter.IntrinsicMatrix(i) == 0
        str = str + "0, ";
    else
        str = str + sprintf("%.5f, ",undistorter.IntrinsicMatrix(i));
    end
end
disp(str + "1]")

disp("Distortion Coeffs");
str = ['[' ...
        num2str( cameraParams.RadialDistortion(1) ) ', '... % k1
        num2str( cameraParams.RadialDistortion(2) ) ', '... % k2
        num2str( cameraParams.TangentialDistortion(1) ) ', '... % p1
        num2str( cameraParams.TangentialDistortion(2) ) ', '... % p1
        num2str( cameraParams.RadialDistortion(3) )... % k3
        ']'];
disp(str);

disp("Optimal matrix:");
str = "[";
for i = 1:8
    if cameraParams.IntrinsicMatrix(i) == 0
        str = str + "0, ";
    else
        str = str + sprintf("%.5f, ",cameraParams.IntrinsicMatrix(i));
    end
end
disp(str + "1]")

%% Save to file
% save('a1001_640_calib', 'undistorter', 'cameraParams');
















%% end
