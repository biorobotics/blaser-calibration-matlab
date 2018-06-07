%% start
% lint three points
P1 = [-9.0807, 7.0378, 20.6196];
P2 = [-40.0592, -31.3176, 87.0532];
P3 = [40.0592, -31.3176, 87.0532];
P4 = [9.0807, 7.0378, 20.6196];
P5 = [1, 1, 1];

%% Solve and Plot
[normal, d] = TriPts2Plane(P1, P3, P4);

disp('ABCD = ');
disp(num2str([normal,d],'\t%.0f'));

figure
x = -100:100; y = -100:100;
[X,Y] = meshgrid(x,y);
Z = (-d - (normal(1)*X) - (normal(2)*Y))/normal(3);
mesh(X,Y,Z)
xlabel('X')
ylabel('Y')
zlabel('Z')
axis equal;
hold on
