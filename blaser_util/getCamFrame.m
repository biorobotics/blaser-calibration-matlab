% X2=7.66; X3= -0.3214 ; Y3= -0.76605;
% 
% 
% X2
% X3
% Y3
% psi = asin(X2/sqrt(1-X3^2))
% theta = asin(-X3)
% phi = asin(Y3/sqrt(1-X3^2))

%% Use Fusion360 to get the axis data

% % Single Blaser 
% % Origin translation vector of the new frame in the old frame
% O_old_in_new = [20.255, 9.783, 5.00]; % [X,Y,Z] in mm
% 
% % tip of new axis vextor project to old axis
% new_axis_length = 20; % mm
% X_tip_old_in_new = [20.255, 9.783, -15.00]; %[X,Y,Z] in mm
% Y_tip_old_in_new = [2.934, -0.217, 5.00]; %[X,Y,Z] in mm
% Z_tip_old_in_new = [10.255, 27.103, 5.00]; %[X,Y,Z] in mm
% 
% X_unit_vec = (X_tip_old_in_new - O_old_in_new) / new_axis_length;
% Y_unit_vec = (Y_tip_old_in_new - O_old_in_new) / new_axis_length;
% Z_unit_vec = (Z_tip_old_in_new - O_old_in_new) / new_axis_length;

% %% Get Euler Angel from Rotation Matrix
% % https://en.wikipedia.org/wiki/Rotation_formalisms_in_three_dimensions#Rotation_matrix_.E2.86.92_Euler_angles_.28z-x-z_extrinsic.29
% 
% R_L = [[-17.321, 7.66, -6.428];
%      [0, 12.856, 15.321];
%      [10,13.268,-11.133]];
% %det(R_L)
%  
% R_R = [[17.321, 7.66, -6.428];
%        [0, -12.856, -15.321];
%        [-10, 13.268, -11.133]]
% %det(R_R)
%    
% R_ref = [[0.5, -0.1465, 0.8535];
%          [0.5, 0.8535, -0.1465];
%          [-0.7072, 0.5, 0.5]];
%    
% R_L = R_L./20;
% R_R = R_R./20;
% % R = R_R./20;
% 
% % (z-x-z extrinsic)
% % phi = atan2(R(3,1),R(3,2))
% % theta = acos(R(3,3))
% % psi = -atan2(R(1,3),R(2,3))
% 
% % Rotation matrix to Quaternion
% % visulasation tool: http://quaternions.online/
% % using Matlab build in func.
% quat_L = rotm2quat(R_L)
% quat_R = rotm2quat(R_R)
% 
% % Wiki Method
% % https://en.wikipedia.org/wiki/Rotation_formalisms_in_three_dimensions#Rotation_matrix_.E2.86.94_quaternion
% % Qr = sqrt(1 + R(1,1) + R(2,2) + R(3,3))/2
% % Qi = (R(3,2) - R(2,3))/(4*Qr)
% % Qj = (R(1,3) - R(3,1))/(4*Qr)
% % Qk = (R(2,1) - R(1,2))/(4*Qr)

%% R Blaser to J8
% Origin translation vector of the new frame in the old frame
Link8_height = 117.564;
Link8_CenterOffset = 25.4;
new_z = Link8_height - Link8_CenterOffset;
O_old_in_new = [33.499, -25.255, 13.867 - new_z]; % [X,Y,Z] in mm
trans_R = O_old_in_new;

% tip of new axis vextor project to old axis
new_axis_length = 20; % mm
X_tip_old_in_new = [20.643, -25.255, -1.454 - new_z]; %[X,Y,Z] in mm
Y_tip_old_in_new = [25.838, -7.934, 20.295 - new_z]; %[X,Y,Z] in mm
Z_tip_old_in_new = [46.767, -15.255, 2.734 - new_z]; %[X,Y,Z] in mm

X_unit_vec = (X_tip_old_in_new - O_old_in_new) / new_axis_length;
Y_unit_vec = (Y_tip_old_in_new - O_old_in_new) / new_axis_length;
Z_unit_vec = (Z_tip_old_in_new - O_old_in_new) / new_axis_length;

% Form Rotation Matrix
R_R = [X_unit_vec; Y_unit_vec; Z_unit_vec]';
quat_R = rotm2quat(R_R);

%% L Blaser to J8
% Origin translation vector of the new frame in the old frame
O_old_in_new = [33.499, 25.255, 13.867 - new_z]; % [X,Y,Z] in mm
trans_L = O_old_in_new;

% tip of new axis vextor project to old axis
new_axis_length = 20; % mm
X_tip_old_in_new = [40.355, 25.255, 29.188 - new_z]; %[X,Y,Z] in mm
Y_tip_old_in_new = [25.838, 7.934, 20.295 - new_z]; %[X,Y,Z] in mm
Z_tip_old_in_new = [46.767, 15.255, 2.734 - new_z]; %[X,Y,Z] in mm

X_unit_vec = (X_tip_old_in_new - O_old_in_new) / new_axis_length;
Y_unit_vec = (Y_tip_old_in_new - O_old_in_new) / new_axis_length;
Z_unit_vec = (Z_tip_old_in_new - O_old_in_new) / new_axis_length;

% Form Rotation Matrix
R_L = [X_unit_vec; Y_unit_vec; Z_unit_vec]';
quat_L = rotm2quat(R_L);

%% Plot rotated frames

%trans_L = [-25.255, 33.499, 13.867];
%trans_R = [25.255, 33.499, 13.867];

figure;
hold on
arrow_length = 10;
oo = [0 0 0];
ex = [1 0 0]*arrow_length;
ey = [0 1 0]*arrow_length;
ez = [0 0 1]*arrow_length;


% plot3([15,0],[0,0],[0,0], 'r'); % X Axis
% plot3([0,0],[15,0],[0,0], 'g'); % Y Axis
% plot3([0,0],[0,0],[15,0], 'b'); % Z Axis

% % axis([-2 2 -2 2 -2 2]);

ux_L = (R_L * ex')';
uy_L = (R_L * ey')';
uz_L = (R_L * ez')';

ux_R = (R_R * ex')';
uy_R = (R_R * ey')';
uz_R = (R_R * ez')';

% ux_test = (R_ref' * ex')';
% uy_test = (R_ref' * ey')';
% uz_test = (R_ref' * ez')';

% Plot reference frame
vectarrow(oo,ex,'r',3);
vectarrow(oo,ey,'g',3);
vectarrow(oo,ez,'b',3);
text(oo(1),oo(2),oo(3),'O_J_8','FontSize',22);

% plot translated and rotated frame L
vectarrow(oo,trans_L,':k',0.5);
vectarrow(trans_L,trans_L+ux_L,'--r',2);
vectarrow(trans_L,trans_L+uy_L,'--g',2);
vectarrow(trans_L,trans_L+uz_L,'--b',2);
text(trans_L(1),trans_L(2),trans_L(3),'O_C_A_M_-_L','FontSize',22);


% plot translated and rotated frame R
vectarrow(oo,trans_R,':k',0.5);
vectarrow(trans_R,trans_R+ux_R,'--r',2);
vectarrow(trans_R,trans_R+uy_R,'--g',2);
vectarrow(trans_R,trans_R+uz_R,'--b',2);
text(trans_R(1),trans_R(2),trans_R(3),'O_C_A_M_-_R','FontSize',22);

% plot rotated frame L
% vectarrow(oo,ux_L,'--r',2);
% vectarrow(oo,uy_L,'--g',2);
% vectarrow(oo,uz_L,'--b',2);

% plot rotated frame R
% vectarrow(oo,ux_R,'--r',2);
% vectarrow(oo,uy_R,'--g',2);
% vectarrow(oo,uz_R,'--b',2);

% plot rotated test frame
% vectarrow(oo,ux_test,'--r',2);
% vectarrow(oo,uy_test,'--g',2);
% vectarrow(oo,uz_test,'--b',2);

%% Display result
str = ['trans_L = [' num2str(trans_L(1)) ', ' num2str(trans_L(2)) ', ' num2str(trans_L(3)) '];' ];
disp(str);

str = ['quat_L = [' num2str(quat_L(1)) ', ' num2str(quat_L(2)) ', ' num2str(quat_L(3)) ', ' num2str(quat_L(4)) '];' ];
disp(str);

str = ['trans_R = [' num2str(trans_R(1)) ', ' num2str(trans_R(2)) ', ' num2str(trans_R(3)) '];' ];
disp(str);

str = ['quat_R = [' num2str(quat_R(1)) ', ' num2str(quat_R(2)) ', ' num2str(quat_R(3)) ', ' num2str(quat_R(4)) '];' ];
disp(str);

disp('matlab quat = [w x y z]');