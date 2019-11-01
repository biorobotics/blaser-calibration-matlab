function [] = show_handeye(A, world, state)
%SHOW_HANDEYE Display handeye calibration results
%   Shows robot arm positions, estimated camera transforms, and the final
%   result location of the checkerboard (which should all line up).
[R, t] = calc_rot_trans(A{1,2}, world, state);
B1 = [R,-R*t; 0,0,0,1];
n_val = size(A,1);

r = state(10);
p = state(11);
w = state(12);
x = state(13);
y = state(14);
z = state(15);

aff = [eul2rotm([w,p,r]),[x;y;z];0,0,0,1];
zz = A{1,1}*aff/B1; % from camera pose to end effector.

figure;
trplot(zz, 'length',.05, 'color','black');
hold on;

n_he = 0;
for i=1:n_val
    if size(A{i,2},1) == 0
        continue
    end
    n_he = n_he+1;
end
cm = hsv(n_he);

% positions = [];
% norms = [];


count = 0;
for i=1:n_val
    if size(A{i,2},1) == 0
        continue
    end
    count = count+1;
    
    [R, t] = calc_rot_trans(A{i,2}, world, state);
    Bi = [R,-R*t; 0,0,0,1];

    % End effector pose
    trplot(A{i,1}, 'length',.04, 'rgb', 'rviz');
    hold on;
    
    % Camera pose (from end effector)
    trplot(A{i,1}*aff, 'length',.05, 'color', cm(count,:));
    hold on;
    
    % Camera pose v2 (from 1st camera)
    camloc = zz*Bi;
    plotCamera('Location', camloc(1:3,4), 'Orientation', camloc(1:3, 1:3)',...
               'Size', .01, 'Color', cm(count,:));
   
    % Checkerboard pose
    trplot(A{i,1}*aff/Bi, 'length',.05, 'color',cm(count,:));
    %disp(A{i,1}*aff/Bi);
   
%    positions = [positions; feat(1:3,4)'];
%    norms = [norms; feat(1:3,3)'];
end
axis tight;
axis vis3d;

end


