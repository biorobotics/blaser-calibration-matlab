addpath('unified_utils');
addpath('blaser_data/640_unified');
clear
close all


n_val = 24;

A = cell(n_val, 3);

figure;
for i = 1:n_val
      fprintf('Image %d\n', i);
      fname = append("im", int2str(i),".png");
      I = imread(fname);
      
      imshow(I);
      hold on;
      drawnow;
    
    laser_pixels = find_laser_new(I);
    %disp(size(laser_pixels));
    
    if numel(laser_pixels) > 1
        lpts = laser_pixels;
        %disp(size(lpts));
        if 1
            disp("Please verify laser stripe integrity");
            plot(lpts(:,1), lpts(:,2), 'g.');
            drawnow;
            
            waitforbuttonpress;

            A{i,3} = lpts;
        end
    end
end
