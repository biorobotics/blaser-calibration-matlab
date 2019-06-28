%% Start
addpath('blaser_util');
addpath('blaser_data');
clear

%%

%camera = 'blaser';
%camera = 'A1003_Calib03';
%camera = 'BlaserA1001_960';
%v = VideoReader([camera,'.mp4']);
%v = VideoReader([camera,'.avi']);

% v = VideoReader('BlaserA1001_960.avi');
v = VideoReader('a1001_640x480_calib_final.mp4');

squareSize = 5; % mm
skip_factor = 9; % integer plz
boardSize = [7,10];
[worldPoints] = generateCheckerboardPoints(boardSize,squareSize);
n_frame = round(v.Duration*v.FrameRate)-2;
disp(n_frame);
threshold = 180;
imagePoints_all = [];
laserPoints_all = [];
frame = [];
% obj = VideoWriter('blaser_calib.avi');
% open(obj)
% for i = 1:n_frame
i = 1;
while i < n_frame
    fprintf('image %d\n', i)
    I = read(v,i);
%     I = read(v,i);
    clf;
    imshow(I);
    drawnow;

    try
        [imagePoints,boardSize_detected] = detectCheckerboardPoints(I);
    catch
        i = i+1;
        continue
    end
    hold on;
    plot(imagePoints(:,1), imagePoints(:,2), 'r');
    drawnow;

    if norm(boardSize - boardSize_detected)>0.1
        i = i+1;
        continue
    end
    try
        imagePoints_all = cat(3,imagePoints_all,flip(imagePoints,1)); %why flip
    catch
        i = i+1;
        continue
    end 
    pixel_data_valid = extractPixelDataFromImg(I, threshold);
    coeffs = polyfit(pixel_data_valid(:,1), pixel_data_valid(:,2), 1);
    [B,fitness] = sort(abs(polyval(coeffs, pixel_data_valid(:,1))-pixel_data_valid(:,2)));
    
    scatter(imagePoints(1,1), imagePoints(1,2), 'p');
    hold on;
    drawnow;

    if (size(fitness, 1) < 100)
        i = i+1;
        continue
    end
    hold on;
    laserPoints_all = cat(3,laserPoints_all,pixel_data_valid(fitness(1:100),:));
    f_end = ceil(size(fitness, 1) * 0.8);
    scatter(pixel_data_valid(fitness(1:f_end),1),pixel_data_valid(fitness(1:f_end),2),5,'g','fill');
    drawnow;
%     scatter(pixel_data_valid(:,1),pixel_data_valid(:,2));
    frame = [frame,i*skip_factor];
%     if i == 1
%         pause();
%     end
%     imwrite(video,[num2str(i),'.tif']);
%     frame = getframe(gcf);
%     writeVideo(obj,frame);
    i = i+skip_factor;

end
% close(obj);
% InitialIntrinsicMatrix = [20000,0,0;0,20000,0;640,360,1];
% [cameraParams,imagesUsed,estimationErrors] = estimateCameraParameters(imagePoints_all,worldPoints,'InitialIntrinsicMatrix',InitialIntrinsicMatrix);
%%
[cameraParams,imagesUsed,estimationErrors] = estimateCameraParameters(imagePoints_all,worldPoints, ...
                                                'NumRadialDistortionCoefficients', 3,  'EstimateTangentialDistortion', true);

%% Check calibration result - LL
figure;
subplot(3,2,[1 3]); showExtrinsics(cameraParams, 'CameraCentric');
subplot(3,2,[2 4]); showExtrinsics(cameraParams, 'PatternCentric');
subplot(3,2,[5 6]);showReprojectionErrors(cameraParams);
displayErrors(estimationErrors, cameraParams);

%% Prepare laser points
K = cameraParams.IntrinsicMatrix';
laser_point = [];
for i = 1:size(frame,2)
    R = cameraParams.RotationMatrices(:,3,i);
    t = cameraParams.TranslationVectors(i,:);
    d = -t*R;
    point_3d = K\[laserPoints_all(:,:,i)';ones(1,size(laserPoints_all,1))];
    z = -d./(point_3d'*R);
    laser_point = [laser_point,point_3d.*[z';z';z']];
end

%% fit laser plane using fit poly11 mechod

% if u see points way out of the norm, kill them here.
laser_point(:, laser_point(2,:) > 0) = [];

[fitresult, gof, plane] = fit3DFitPlane(laser_point(1,:),laser_point(2,:),laser_point(3,:));

%% Prepare laser points and fit plane using SVD method
% K = cameraParams.IntrinsicMatrix';
% laser_point = [];
% for i = 1:size(frame,2)
%     R = cameraParams.RotationMatrices(:,3,i);
%     t = cameraParams.TranslationVectors(i,:);
%     d = -t*R;
%     point_3d = K\[laserPoints_all(:,:,i)';ones(1,size(laserPoints_all,1))];
%     z = -d./(point_3d'*R);
%     laser_point = [laser_point,point_3d.*[z';z';z']];
% end
% 
% % homogenous form of laser points
% pts = [laser_point', ones(size(laser_point,2),1)];
% 
% % homogenous form of camera intrinsic matrix
% % from K = 3*3 to K = 4*4
% K_expand = [K,zeros(3,1);zeros(1,3),1];
% AM = pts * (K_expand ^ (-1))';
% [U,S,V] = svd(AM);
% 
% gc = V(:,4)';
% 
% % gc is the camera parameters needed from distance detection
% % it can be verified that gc(1)^2+gc(2)^2+gc(3)^2 = 1
% gc_square_sum = sqrt(gc(1)^2+gc(2)^2+gc(3)^2);
% gc_gain = 1/gc_square_sum;
% gc = V(:,4) * gc_gain;
% gc = gc'
% 
% %% plot calulated laser plane
% hold on
% A = gc(1); B = gc(2); C = gc(3); D = gc(4);
% point = [0,0,-D/C];
% normal = [A,B,C];
% 
% %# a plane is a*x+b*y+c*z+d=0
% %# [a,b,c] is the normal. Thus, we have to calculate
% %# d and we're set
% d = -point*normal'; %'# dot product for less typing
% 
% %# create x,y
% [xx,yy]=ndgrid(-40:40,-15:0);
% 
% %# calculate corresponding z
% z = (-normal(1)*xx - normal(2)*yy - d)/normal(3);
% 
% %# plot the surface
% surf(xx,yy,z)
% hold off

%% plot laser points
figure
hold on 
scatter3(laser_point(1,:),laser_point(2,:),laser_point(3,:));

%% show extrinsics and laser points
figure
showExtrinsics(cameraParams, 'CameraCentric');
hold on 
scatter3(laser_point(1,:),laser_point(3,:),laser_point(2,:),'r');
%plot( fitresult );
%plot( fitresult, [plane.xData, plane.zData], plane.yData );

%% polot a surface
%# create x,y
[xx,yy] = meshgrid(-20:20,-20:20);

%# calculate corresponding z
zz = -1*(plane.a*xx + plane.b*yy + plane.d)/plane.c;

text = ['plane.a = ' num2str(plane.a)];  disp(text); 
text = ['plane.b = ' num2str(plane.b)];  disp(text); 
text = ['plane.c = ' num2str(plane.c)];  disp(text); 
text = ['plane.d = ' num2str(plane.d)];  disp(text); 

text = ['laser_plane (allx100) = [' num2str(plane.a*100) ', ' num2str(plane.b*100) ', ' num2str(plane.c*100) ', ' num2str(plane.d*100) ']'];
disp(text); 

%# plot the surface
hold on
surf(xx,zz,yy)
hold off

%% polot a ideal surface when no z rotation

figure
showExtrinsics(cameraParams, 'CameraCentric');
hold on 
scatter3(laser_point(1,:),laser_point(3,:),laser_point(2,:),'r');
%plot( fitresult );
%plot( fitresult, [plane.xData, plane.zData], plane.yData );

%# create x,y
[xx,yy] = meshgrid(-20:20,-20:20);

%# calculate corresponding z
zz = -1*(plane.a/5*xx + plane.b*yy + plane.d)/plane.c;

text = ['plane.a = ' num2str(plane.a/5)];  disp(text); 
text = ['plane.b = ' num2str(plane.b)];  disp(text); 
text = ['plane.c = ' num2str(plane.c)];  disp(text); 
text = ['plane.d = ' num2str(plane.d)];  disp(text); 

%# plot the surface
hold on
surf(xx,zz,yy)
hold off

%% Print result
disp('------------------------------------------------------------------'); 
%cameraParams.PrincipalPoint
% camera matrix
% [fx 0 cx, 0 fy cy, 0 0 1]
text = ['camera matrix = [' ...
    num2str( cameraParams.IntrinsicMatrix(1,1) ) ', '... % fx
    num2str( cameraParams.IntrinsicMatrix(2,1) ) ', '... % 0
    num2str( cameraParams.IntrinsicMatrix(3,1) ) ', '... % cx
    num2str( cameraParams.IntrinsicMatrix(1,2) ) ', '... % fy
    num2str( cameraParams.IntrinsicMatrix(2,2) ) ', '... % 0
    num2str( cameraParams.IntrinsicMatrix(3,2) ) ', '... % cy
    '0, 0, 1]'];
disp('camera matrix = [fx 0 cx, 0 fy cy, 0 0 1]');
disp(text);

% distortion_coefficients
text = ['distortion_coefficients = [' ...
        num2str( cameraParams.RadialDistortion(1) ) ', '... % k1
        num2str( cameraParams.RadialDistortion(2) ) ', '... % k2
        num2str( cameraParams.TangentialDistortion(1) ) ', '... % p1
        num2str( cameraParams.TangentialDistortion(2) ) ', '... % p1
        num2str( cameraParams.RadialDistortion(3) )... % k3
        ']'];
disp('distortion_coefficients = [k1 k2 p1 p2 k3] - k:radial, p:tangential');
disp(text);

% Laser Plane 
text = ['laser_plane (allx100) = [' num2str(plane.a*100) ', ' num2str(plane.b*100) ', ' num2str(plane.c*100) ', ' num2str(plane.d*100) ']'];
disp(text); 

disp('------------------------------------------------------------------');

%% Compare with 'ground truth'

figure
% showExtrinsics(cameraParams, 'CameraCentric');
scatter3(laser_point(1,:),laser_point(3,:),laser_point(2,:),'r');

[xx,yy] = meshgrid(-20:20,-30:20);
zz = -1*(plane.a/5*xx + plane.b*yy + plane.d)/plane.c;
hold on
surf(xx,zz,yy);

zz_gt = -1*(6.4656*xx -209.0053*yy + 3554.7677)/(-100);
hold on
s = surf(xx,zz_gt,yy);
s.CData = 20+ones(size(yy));

% figure
% surf(xx, yy, zz-zz_gt)

%% Hand Eye Calibration

%img_hand = imread('BlaserEndEffector.png');
img_hand = imread('blaser_1280x960.png');

[hand_imagePoints,hand_boardSize] = detectCheckerboardPoints(img_hand);

hand_squareSize = 6; % in millimeters
hand_worldPoints = generateCheckerboardPoints(hand_boardSize, hand_squareSize);

% show original image
figure;
subplot(1,2,1);
imshow(img_hand);
hold on;
plot(hand_imagePoints(:,1),hand_imagePoints(:,2),'ro');

% undistort image
[im, newOrigin] = undistortImage(img_hand, cameraParams, 'OutputView', 'full');
subplot(1,2,2); imshow(im);
hold on;
%plot(hand_imagePoints(:,1),hand_imagePoints(:,2),'ro');
[hand_imagePoints,hand_boardSize] = detectCheckerboardPoints(im);

%
hand_imagePoints = [hand_imagePoints(:,1) + newOrigin(1), ...
             hand_imagePoints(:,2) + newOrigin(2)];

plot(hand_imagePoints(:,1),hand_imagePoints(:,2),'bo');
plot(hand_imagePoints(1,1),hand_imagePoints(1,2),'r+');

[hand_R, hand_t] = extrinsics(hand_imagePoints, hand_worldPoints, cameraParams);

% Compute camera pose.
[orientation, location] = extrinsicsToCameraPose(hand_R, hand_t);

figure
plotCamera('Location',location,'Orientation',orientation,'Size',5);
hold on
pcshow([worldPoints,zeros(size(worldPoints,1),1)], ...
  'VerticalAxisDir','down','MarkerSize',40);

imagePoints_origin = hand_imagePoints(1,:);
worldPoints_origin = pointsToWorld(cameraParams, hand_R, hand_t, imagePoints_origin);

worldPoints_origin = [worldPoints_origin 0];
[~, cameraLocation] = extrinsicsToCameraPose(hand_R, hand_t);
distanceToCamera = norm(worldPoints_origin - cameraLocation);

disp(hand_R)
disp(hand_t)

%% Save to file
% fileID = fopen([camera,'.yaml'],'w');
% fprintf(fileID,'Camera intrisic matrix:\n');
% fprintf(fileID,'%f %f %f\n',cameraParams.IntrinsicMatrix);
% fprintf(fileID,'laser plane in camera frame:\n');
% %fprintf(fileID,'%f %f %f %f\n',gc);
% fclose(fileID);
