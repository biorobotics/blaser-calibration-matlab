% function [retpts, jac_x, jac_y] = undistort_points(pts, im_sz, cx, cy, k1, k2, p1, p2, k3)
function [retpts] = undistort_points(pts, state)
% Undistorts pts according to Brown-Conrady method.
% pts a Nx2 vector
% all the others scalars
% jac_x returns Nx7 (cx, cy, k1, k2, p1, p2, k3)
% jac_y returns Nx7 (cx, cy, k1, k2, p1, p2, k3)

    fx = state(1);
    fy = state(2);
    cx = state(3);
    cy = state(4);
    k1 = state(5);
    k2 = state(6);
    p1 = state(7);
    p2 = state(8);
%     p1 = 0; p2 = 0;
    k3 = state(9);

    K = [fx,0,0;0,fy,0;cx,cy,1];
    cp = cameraParameters('IntrinsicMatrix',K,...
        'RadialDistortion',[k1,k2,k3],'TangentialDistortion',[p1,p2]);
    
    % Inverse radial distortion coeffs
    % source: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4934233/
    % doi: 10.3390/s16060807
%     b1 = -k1;
%     b2 = 3*k1^2 - k2;
%     b3 = -12*k1^3 + 8*k1*k2 - k3;
%     b4 = 55*k1^4 - 55*k2*k1^2 + 5*k2^2 + 10*k1*k3;
%     
%     x0 = pts(:,1) - cx;
%     y0 = pts(:,2) - cy;
%     rsq = x0.^2 + y0.^2;
%     
%     cdist = rsq.*(b1 + rsq.*(b2 + rsq.*(b3 + rsq.*b4)));
%     xu = pts(:,1) + x0.*cdist;
%     yu = pts(:,2) + y0.*cdist;

    retpts = cp.undistortPointsImpl(pts);

%     dx = (0:360)';
%     inp = [cx+dx,cy*ones(size(dx))];
%     blyats = cp.undistortPointsImpl(inp);
%     
% %     retpts = undistortPoints(pts, cp);
%     retpts = cp.undistortPointsImpl(pts);
%     retpts2 = [xu,yu];
%     diff = retpts - retpts2;
    
%     x0 = (pts(:,1) - cx) / fx;
%     x0 = pts(:,1);
% %     y0 = (pts(:,2) - cy) / fy;
%     y0 = pts(:,2);
%     rsq = x0.^2 + y0.^2;
%     
%     % formula according to Brown-Conrady
% %     rad_distorts = k1 * rsq + k2*rsq.^2 + k3*rsq.^3;
%     cdist = 1 + ((k3*rsq + k2).*rsq + k1).*rsq;
%     
%     dx = p2 * (rsq + 2*x0.^2) + 2*p1 * x0.*y0;
%     dy = p1 * (rsq + 2*y0.^2) + 2*p2 * x0.*y0;
%     
%     xu = x0 .* cdist + dx;
%     yu = y0 .* cdist + dy;
% %     xu = pts(:,1) + x0./cdist + dx;
% %     yu = pts(:,2) + y0./cdist + dy;
%         
%     
%     retpts = [xu, yu];
    
%     % jac of radial distortions
%     drdcx = -2*k1*x0 - 4*k2*x0.*r_sq - 6*k3* ((r_sq.^2) .* x0);
%     drdcy = -2*k1*y0 - 4*k2*y0.*r_sq - 6*k3* ((r_sq.^2) .* y0);
%     drdk1 = r_sq;
%     drdk2 = r_sq.^2;
%     drdk3 = r_sq.^3;
%     
%     jxcx = x0.*drdcx - rad_distorts -6*p1*x0 - p2*y0;
%     jxcy = x0.*drdcy - 2*p1*y0 - p2*x0;
%     jac_x = [jxcx/w, ...
%              jxcy/h, ...
%              x0.*drdk1, ...
%              x0.*drdk2, ...
%              (r_sq + 2*x0.^2), ...
%              x0.*y0, ...
%              x0.*drdk3];
%          
%     jycx = y0.*drdcx - 2*p2*x0 - p1*y0;
%     jycy = y0.*drdcy - rad_distorts -6*p2*y0 - p1*x0;
%     jac_y = [jycx/w, ...
%              jycy/h, ...
%              y0.*drdk1, ...
%              y0.*drdk2, ...
%              x0.*y0, ...
%              (r_sq + 2*y0.^2), ...
%              y0.*drdk3];
end
