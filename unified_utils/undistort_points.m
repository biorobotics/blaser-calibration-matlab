function [retpts, Jrx, Jry] = undistort_points(pts, state)
% Undistorts pts according to Brown-Conrady method.
% pts a Nx2 vector
% all the others scalars
% Jrx returns Nx9 (fx, fy, cx, cy, k1, k2, p1, p2, k3)
% Jry returns Nx9 (fx, fy, cx, cy, k1, k2, p1, p2, k3)

    fx = state(1);
    fy = state(2);
    cx = state(3);
    cy = state(4);
    k1 = state(5);
    k2 = state(6);
    p1 = state(7);
    p2 = state(8);
    k3 = state(9);
    npt = size(pts,1);
    
%     syms fx fy cx cy k1 k2 p1 p2 k3;
% cp = cameraParameters('IntrinsicMatrix',K,'RadialDistortion',[k1,k2,k3],'TangentialDistortion',[p1,p2]);
% cp = cameraParameters('IntrinsicMatrix',eye(3),'RadialDistortion',[k1,k2,k3],'TangentialDistortion',[p1,p2]);
% cp = cameraParameters('IntrinsicMatrix',eye(3),'RadialDistortion',[k1,k2,k3]);
    
    K = [fx,0,0;0,fy,0;cx,cy,1];
    Ki = [1/fx,0,0;0,1/fy,0;-cx/fx,-cy/fy,1];

    % forward distortion
    zpts = [pts,ones(npt,1)] * Ki;
    x0 = zpts(:,1);
    y0 = zpts(:,2);
    rsq = x0.^2 + y0.^2;
    cdist = 1 - rsq.*(k1 + rsq.*(k2 + k3*rsq));
    dx = p2 * (rsq + 2*x0.^2) + 2*p1 * x0.*y0;
    dy = p1 * (rsq + 2*y0.^2) + 2*p2 * x0.*y0;

    xg = x0.*cdist - dx;
    yg = y0.*cdist - dy;
    rsqg = xg.^2 + yg.^2;
    rad_g = 1 + rsqg.*(k1 + rsqg.*(k2 + k3*rsqg));
    tanx_g = p2 * (rsqg + 2*xg.^2) + 2*p1 * xg.*yg;
    tany_g = p1 * (rsqg + 2*yg.^2) + 2*p2 * xg.*yg;
    xres = xg.*rad_g + tanx_g - x0;
    yres = yg.*rad_g + tany_g - y0;
%     resid = [xres,yres];

    % Jacobian of the error function (distortion)
    Jxx = 1 + 6*p2*xg + p1*yg + k1*(2*xg.^2+rsqg) + rsqg.*(k2*(4*xg.^2+rsqg) + k3*rsqg.*(6*xg.^2+rsqg));
    Jxy = p1*xg + 2*yg.*(p2 + xg.*(k1 + rsqg.*(2*k2 + 3*k3*rsqg)));
    Jyx = p2*yg + 2*xg.*(p1 + yg.*(k1 + rsqg.*(2*k2 + 3*k3*rsqg)));
    Jyy = 1 + p2*xg + 6*p1*yg + k1*(2*yg.^2+rsqg) + rsqg.*(k2*(4*yg.^2+rsqg) + k3*rsqg.*(6*yg.^2+rsqg));
    
    % Get jacobian inverse
    det = Jxx.*Jyy - Jxy.*Jyx;
%     Jix = [Jyy, -Jxy];
%     Jiy = [-Jyx, Jxx];
    
    % Solve?
    xu = xg - (Jyy.*xres - Jxy.*yres)./det;
    yu = yg - (Jxx.*yres - Jyx.*xres)./det;
    redist = [xu, yu, ones(npt,1)]*K;
    retpts = redist(:,1:2);
    return

%     xg = x0.*cdist;
%     yg = y0.*cdist;
%     rg = r.*cdist;
    rsqg = rsq.*cdist.^2;
    % Taylor series expansion 1st order
    c0 = cdist.*(1+rsqg.*(k1 + rsqg.*(k2 + k3*rsqg))) - 1;
    c1 = 1 + rsqg.*(3*k1 + rsqg.*(5*k2 + 7*k3*rsqg));
    rad_undist = cdist - c0./c1;


    xd = x0.*rad_undist - dx;
    yd = y0.*rad_undist - dy;
    redist = [xd, yd, ones(npt,1)]*K;
    retpts = redist(:,1:2);
%     retpts = redist(:,1);

    if nargout == 1
        return
    end

    % Jacobian
    Jx = [-x0/fx, zeros(size(x0)),  -ones(size(x0))/fx, zeros(size(x0))];
    Jy = [zeros(size(y0)), -y0/fy, zeros(size(x0)), -ones(size(y0))/fy ];
    Jrsq = 2*(x0.*Jx + y0.*Jy);

    Jcdist = -(k1 + rsq.*(2*k2 + 3*k3*rsq)).*Jrsq;
    Jk1 = -rsq;
    Jk2 = -rsq.^2;
    Jk3 = -rsq.^3;

    Jdx = p2*(Jrsq + x0.*Jx*4) + 2*p1*(x0.*Jy + y0.*Jx);
    Jdy = p1*(Jrsq + y0.*Jy*4) + 2*p2*(x0.*Jy + y0.*Jx);

    Jxd = [x0.*Jcdist + cdist.*Jx - Jdx,...
           x0.*Jk1, x0.*Jk2, -2*x0.*y0, -rsq - 2*x0.^2, x0.*Jk3];

    Jyd = [y0.*Jcdist + cdist.*Jy - Jdy,...
           y0.*Jk1, y0.*Jk2, -rsq - 2*y0.^2, -2*x0.*y0, y0.*Jk3];

    Jrx = fx*Jxd;
    Jrx(:,3) = Jrx(:,3)+1;
    Jrx(:,1) = Jrx(:,1) + xd;
    
    Jry = fy*Jyd;
    Jry(:,4) = Jry(:,4)+1;
    Jry(:,2) = Jry(:,2) + yd;


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


    % inverse radial coefficients
%     b1 = -k1;
%     b2 = 3*k1^2 - k2;
%     b3 = -12*k1^3 + 8*k1*k2 - k3;
%     b4 = 55*k1^4 - 55*k2*k1^2 + 5*k2^2 + 10*k1*k3;
%     b5 = -273*k1^5 + 364*k2*k1^3 - 78*k1*k2^2 - 78*k3*k1^2 + 12*k2*k3;
%     b6 = 1428*k1^6 - 2380*k1^4*k2 + 840*k1^2*k2^2 - 35*k2^3 + ...
%         560*k1^3*k3 - 210*k1*k2*k3 + 7*k3^2;
%     b7 = -7752*k1^7 + 15504*k1^5*k2 - 7752*k1^3*k2^2 + 816*k1*k2^3 - ...
%         3876*k1^4*k3 + 2448*k1^2*k2*k3 - 136*k2^2*k3 - 136*k1*k3^2;

end
