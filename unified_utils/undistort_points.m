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
    npt = size(pts,1);

    K = [fx,0,0;0,fy,0;cx,cy,1];
%     cp = cameraParameters('IntrinsicMatrix',K,...
%         'RadialDistortion',[k1,k2,k3],'TangentialDistortion',[p1,p2]);
%     redist = cp.undistortPointsImpl(pts);
%     redist = redist(:,1:2);
    

    % forward distortion
    zpts = [pts,ones(npt,1)] / K;
    x0 = zpts(:,1);
    y0 = zpts(:,2);
    rsq = x0.^2 + y0.^2;
    cdist = 1 + rsq.*(k1 + rsq.*(k2 + k3*rsq));
    
    dx = p2 * (rsq + 2*x0.^2) + 2*p1 * x0.*y0;
    dy = p1 * (rsq + 2*y0.^2) + 2*p2 * x0.*y0;
    
    xd = x0.*cdist + dx;
    yd = y0.*cdist + dy;
    redist = [2*x0-xd, 2*y0-yd, ones(npt,1)]*K;
    retpts = redist(:,1:2);

    % Inverse radial distortion coeffs
    % source: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4934233/
    % doi: 10.3390/s16060807
%     b1 = -k1;
%     b2 = 3*k1^2 - k2;
%     b3 = -12*k1^3 + 8*k1*k2 - k3;
%     b4 = 55*k1^4 - 55*k2*k1^2 + 5*k2^2 + 10*k1*k3;
%     b5 = -273*k1^5 + 364*k2*k1^3 - 78*k1*k2^2 - 78*k3*k1^2 + 12*k2*k3;
%     b6 = 1428*k1^6 - 2380*k1^4*k2 + 840*k1^2*k2^2 - 35*k2^3 + ...
%          560*k1^3*k3 - 210*k1*k2*k3 + 7*k3^2;
%     
%     zpts = [pts,ones(npt,1)] / K;
%     x0 = zpts(:,1);
%     y0 = zpts(:,2);
%     rsq = x0.^2 + y0.^2;
%     
%     inv_cdist = 1 + rsq.*(b1 + rsq.*(b2 + rsq.*(b3 + rsq.*(b4 + rsq.*(b5 + rsq.*b6)))));
%     xu_r = x0.*inv_cdist;
%     yu_r = y0.*inv_cdist;
    
    % Undistort happened assuming tangential is negligible.
%     rsq_r = xu_r.^2 + yu_r.^2;
%     dx_t = p2 * (rsq_r + 2*xu_r.^2) + 2*p1 * xu_r.*yu_r;
%     dy_t = p1 * (rsq_r + 2*yu_r.^2) + 2*p2 * xu_r.*yu_r;
    
    % Since tangential small, assume it's locally linear (same forward and back)
%     x0 = zpts(:,1) - dx_t;
%     y0 = zpts(:,2) - dy_t;
%     rsq = x0.^2 + y0.^2;
%     inv_cdist = 1 + rsq.*(b1 + rsq.*(b2 + rsq.*(b3 + rsq.*(b4 + b5*rsq))));
%     xu = x0.*inv_cdist;
%     yu = y0.*inv_cdist;
%     undistorted = [xu_r-dx_t,yu_r-dy_t,ones(npt,1)]*K;
%     retpts = undistorted(:,1:2);
    
%     x0 = pts(:,1) - cx;
%     y0 = pts(:,2) - cy;
%     rsq = x0.^2 + y0.^2;
%     
%     cdist = rsq.*(b1 + rsq.*(b2 + rsq.*(b3 + rsq.*b4)));
%     xu = pts(:,1) + x0.*cdist;
%     yu = pts(:,2) + y0.*cdist;

end
