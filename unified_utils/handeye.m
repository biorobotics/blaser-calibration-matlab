function [resid, ch_loc, ch_ori] = handeye(A, world, state)
%HANDEYE Does hand-eye calibration via AX=XB
%   Include reprojection in here to speed up.
% A- cell of inputs (robot arm pose + camera points)
% world- world points (used for reprojection errors lmao)
% state- 1x15 vector [fx fy cx cy k1 k2 p1 p2 k3 r p w x y z]

n_im = size(A,1);

r = state(10);
p = state(11);
w = state(12);
x = state(13);
y = state(14);
z = state(15);
rot = eul2rotm([w,p,r]);
aff = [rot,[x;y;z];0,0,0,1];

positions = [];
norms = [];

count = 0;
AB = cell(n_im, 2);

resid = [];
for i=1:n_im
    if size(A{i,2}, 1) == 0
        continue
    end
    
    count = count+1;
    
    [err, R, t] = reproj_err(A{i,2}, world, state);
    resid = [resid; err/10];
    
    B = [R,-R*t; 0,0,0,1];
    
    AB{count, 1} = A{i,1};
    AB{count, 2} = B;
    
    feat = A{i,1}*aff/B;
    positions = [positions; feat(1:3,4)'];
    norms = [norms; feat(1:3,3)'];
end

ch_loc = mean(positions)';
ch_ori = mean(norms)';

handeye_resid = [];
for i=1:count
    for j=1:count
        if i==j
            continue
        end
        da = AB{i,1}\AB{j,1};
        db = AB{i,2}\AB{j,2};
        
        Ra = da(1:3,1:3);
        Rb = db(1:3,1:3);
        ta = da(1:3,4);
        tb = db(1:3,4);
        
        F = Ra * rot - rot * Rb;
        F = F(:);
        G = (Ra - eye(3))*[x;y;z] - rot*tb + ta;
        
        handeye_resid = [handeye_resid;F;G];
    end
end

resid = [resid; handeye_resid];

end

