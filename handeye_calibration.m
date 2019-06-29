%% Start
addpath('blaser_data/handeye_data');
clear

hand_boardSize = [7 10];
hand_squareSize = .005; %metess
hand_worldPoints = generateCheckerboardPoints(hand_boardSize, hand_squareSize);

%% load file
f = fopen('data.txt','r');
n_val = 9;
A = cell(n_val, 1);
B = cell(n_val, 1);

imagePoints_all = [];
laserPoints_all = [];

for i = 1:n_val
    tr = fscanf(f, '[%f, %f, %f]\n', 3);
    rot = fscanf(f, '[%f, %f, %f, %f]\n', 4);
    quat = [rot(4), rot(1), rot(2), rot(3)];
    A{i} = [quat2rotm(quat), tr; 0,0,0,1];
    
    fname = fscanf(f, '%s\n', 1);
    img = imread(fname);
    [imagePoints,boardSize_detected] = detectCheckerboardPoints(img);
    imagePoints_all = cat(3,imagePoints_all,flip(imagePoints,1)); %why flip
end

%% do camera calibration from images
[cameraParams,imagesUsed,estimationErrors] = estimateCameraParameters(imagePoints_all,hand_worldPoints, 'worldUnits', 'm', ...
                                                'NumRadialDistortionCoefficients', 3,  'EstimateTangentialDistortion', true);
for i=1:n_val
    cam_rot = cameraParams.RotationMatrices(:,:,i);
    cam_tr  = cameraParams.TranslationVectors(i,:);
    
    B{i} = [cam_rot,-cam_rot*cam_tr'; 0,0,0,1];
end

%% visualize loaded thingys
% Tool pose (x down)
figure;
for i = 1:n_val
    trplot(A{i}, 'length',.1, 'rgb', 'rviz');
    hold on;
end
axis tight;
axis vis3d;

% camera pose (z down)
figure
showExtrinsics(cameraParams, 'PatternCentric');
hold on;
for i=1:n_val
    trplot(B{i}, 'length', .02, 'rgb');
end

axis tight;
axis vis3d;

%% Optimization-based calibration method (init step)
% https://pdfs.semanticscholar.org/d36c/d02f387724687ac3276f4cc449a30908c710.pdf
% iterative optimization formula
da = A{1}\A{2};
db = B{1}\B{2};
sing_eq = kron(da(1:3,1:3), eye(3)) - kron(eye(3), db(1:3,1:3)');
[ev,D] = eigs(sing_eq'*sing_eq, 1, 'sm');
[U,S,V] = svd(vec2mat(ev,3));
R = U*V';
R = R*sign(det(R));
% t = (da(1:3,1:3)-eye(3))\(R*db(1:3,4)-da(1:3,4));
X = [rotm2eul(R)';0;0;0]; % the math for initial position guess wasn't working out

%% iter step setup
syms x y z r p w
assume(x, 'real')
assume(y, 'real')
assume(z, 'real')
assume(r, 'real')
assume(p, 'real')
assume(y, 'real')

rot = [cos(w),-sin(w),0;sin(w),cos(w),0;0,0,1]...
    *[cos(p),0,sin(p);0,1,0;-sin(p),0,cos(p)]...
    *[1,0,0;0,cos(r),-sin(r);0,sin(r),cos(r)];

H = [];
J = [];
f_w = 1; %positional weight
g_w = 1; %angular weight
for i=1:n_val
    for j=1:n_val
        if i==j
            continue
        end
        da = A{i}\A{j};
        db = B{i}\B{j};
        
        Ra = da(1:3,1:3);
        Rb = db(1:3,1:3);
        ta = da(1:3,4);
        tb = db(1:3,4);
        
        F = Ra * rot - rot * Rb;
        F = f_w * F(:);
        G = g_w * (Ra - eye(3))*[x;y;z] - rot*tb + ta;
        
        Ji = jacobian([F;G], [r,p,y,x,y,z]);

        H = [H;F;G];
        J = [J;Ji];
    end
end

jtj = combine(J' * J);
jth = combine(J' * H);

%% iter step (runs until error decreases by less than 1E-7)
prev_err = inf;
cll = num2cell(X);
[r,p,w,x,y,z] = cll{:};
err = norm(vpa(subs(H)));
disp('Init error')
disp(err)

% Using levenberg marquardt
mu = .8;
while prev_err/err-1 > 1E-7
    lhs = vpa(subs(jtj)) + mu*eye(6);
    rhs = -vpa(subs(jth));
    dX = double(lhs\rhs);
    X_prev = X;
    X = X+dX;

    prev_err = err;
    cll = num2cell(X);
    [r,p,w,x,y,z] = cll{:};
    err = norm(vpa(subs(H)));
    disp('Error:');
    disp(err);
end
disp('Done!');
disp('Final err:');
disp(err);
beep;

%% visualize the optimized transform
aff = [eul2rotm([X(3) X(2) X(1)]),X(4:6);0,0,0,1];
figure;
zz = A{1}*aff/(B{1});

trplot(zz, 'length',.05, 'color','black');
hold on;
cm = hsv(n_val);
err = [];
for i=1:n_val
    trplot(A{i}, 'length',.04, 'rgb', 'rviz');
    hold on;
    trplot(A{i}*aff, 'length',.05, 'color', cm(i,:));
    hold on;
    
    camloc = zz*B{i};
    plotCamera('Location', camloc(1:3,4), 'Orientation', camloc(1:3, 1:3)',...
               'Size', .01, 'Color', cm(i,:));
   
   feat = A{i}*aff/B{i};
   trplot(A{i}*aff/B{i}, 'length',.05, 'color',cm(i,:));
   err = [err; feat(1:3,4)'];
end
axis tight;
axis vis3d;

disp(std(err))
% print in format easy for URDF
sprintf('xyz="%f %f %f" rpy="%f %f %f"', X(4), X(5), X(6), X(1), X(2), X(3))
