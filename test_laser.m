addpath('unified_utils');
addpath('blaser_data/640_unified');
%addpath('/home/haowensh/Projects/blaser/calibration/laser_embedded')
clear
close all


n_val = 24;
minLaserStripePts = 70;

A = cell(n_val, 3);

figure;
for ii = 1:n_val
      fprintf('Image %d\n', ii);
      fname = append("im", int2str(ii),".png");
      I = imread(fname);
      
      imshow(I);
      hold on;
      drawnow;
    
    laser_pixels = find_laser_new(I);
    %disp(size(laser_pixels));
    
    if numel(laser_pixels) > 1
        coeffs = polyfit(laser_pixels(:,1), laser_pixels(:,2), 1);
        dists = abs(polyval(coeffs, laser_pixels(:,1)) - laser_pixels(:,2));
        lpts = laser_pixels(dists < 2, :);

        coeffs2 = polyfit(lpts(:,1), lpts(:,2), 1);
        dists = abs(polyval(coeffs2, laser_pixels(:,1)) - laser_pixels(:,2));
        lpts = laser_pixels(dists < 1, :);

        if size(lpts, 1) > minLaserStripePts
            disp("Laser line found!");
            hold on
            % Plot found points
            plot(lpts(:,1), lpts(:,2), 'g.');
            % Plot interpolated line
            X = linspace(1, size(I, 2), size(I, 2));
            fline = polyval(coeffs, X);
            plot(X, fline, '--', 'color', [0 1 1])
            hold off
            drawnow;
            waitforbuttonpress;

            A{ii,3} = lpts;
        end
    end
end
