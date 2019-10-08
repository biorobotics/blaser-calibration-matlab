% img  = imread('blaser_data/1280_unified/im0.png');
% find_laser(img);


function [pts] = find_laser_new(img)
    pts = [];
%     img_copy = img;
%     figure;
%     imshow(img);
%     r = img(:,:,1);
%     empty_b = zeros(size(r));
%     empty_g = zeros(size(r)); 
%     only_red = cat(3,r, empty_g, empty_b);
%     figure;

%     Median Blur
    img = medfilt3(img);
    ksize = [7 7];
    sigma = [0.3*((ksize(2)-1)*0.5 - 1) + 0.8, 0.3*((ksize(1)-1)*0.5 - 1) + 0.8];
%     Gaussian Blur
    img = imgaussfilt(img,sigma, 'FilterSize',ksize);
    gb1 = int32(img(:,:,3));
    
%     Color Balance
    im_correct = color_balance(img,0.0, 0);
%     imshow(im_correct);
%       Masks    
    rmask = getredmask(im_correct);
    smask = getsatmask(im_correct);
    tmask = rmask+smask;
   
    im_comb = bitwise_and(gb1,tmask);
    m_img = mean(im_comb, 'all');
    im_sub_mean = im_comb - m_img;
    im_sub_mean(im_sub_mean <0)=0;
%     figure;
%     imshow(im_sub_mean);
    
%     max_pixel = max(im_sub_mean, [], 'all');
    [val_max, col_max] = max(im_sub_mean, [], 1); 
    for i=1:size(img,1)
        j = col_max(i)-1;
   
        while (j >= 1 && im_sub_mean(j,i) > 0.8*val_max(i))
            j=j-1;
        end
        k = col_max(i)+1;
        
        while (k <= size(img,1) && im_sub_mean(k,i) > 0.8*val_max(i))
            k=k+1;
        end
        weighted_sum = 0.0;
        val_sum = 0.0;
        for idx=j+1:k-1
            weighted_sum = weighted_sum+ im_sub_mean(idx,i)*(idx-j);
            val_sum = val_sum +im_sub_mean(idx,i);
        end
        if (val_sum <=0)
            continue;
        end
        pts = [pts; i, double(weighted_sum/val_sum + j)];
    end
%     figure; 
%     imshow(img_copy);
%     hold on;
%     drawnow; 
%     plot(pts(:,1), pts(:,2), 'g.');
%     drawnow;
end


function [res] = bitwise_and(img,mask)
    res = zeros(size(img)); 
    for i=1:size(img,1)
        for j=1:size(img,2)
            if (mask(i,j) ~= 0)
                res(i,j)= bitand(img(i,j), img(i,j));
            end
        end
    end

end


function [mask] = getredmask(img)

    % Convert RGB image to chosen color space
    hsv_img = rgb2hsv(img);
    % Define thresholds for channel 1 based on histogram settings
    channel1Min = 0.796;
    channel1Max = 0.166;
    % Define thresholds for channel 2 based on histogram settings
    channel2Min = 0.001;
    channel2Max = 1.00;
    % Define thresholds for channel 3 based on histogram settings
    channel3Min = 0.588;
    channel3Max = 1.000;
    % Create mask based on chosen histogram thresholds
    sliderBW = ( (hsv_img(:,:,1) >= channel1Min) | (hsv_img(:,:,1) <= channel1Max) ) & ...
        (hsv_img(:,:,2) >= channel2Min ) & (hsv_img(:,:,2) <= channel2Max) & ...
        (hsv_img(:,:,3) >= channel3Min ) & (hsv_img(:,:,3) <= channel3Max);
    mask = sliderBW;

    mask = imdilate(mask, strel('rectangle',[9 9])); 
    mask = imclose(mask, strel('rectangle',[25 35]));

    
end

function [mask] = getsatmask(img)

    hsv_img = rgb2hsv(img);
    channel1Min = 0.000;
    channel1Max = 1.000;

    % Define thresholds for channel 2 based on histogram settings
    channel2Min = 0.000;
    channel2Max = 0.227;

    % Define thresholds for channel 3 based on histogram settings
    channel3Min = 0.827;
    channel3Max = 1.000;

    % Create mask based on chosen histogram thresholds
    sliderBW = ( (hsv_img(:,:,1) >= channel1Min) | (hsv_img(:,:,1) <= channel1Max) ) & ...
        (hsv_img(:,:,2) >= channel2Min ) & (hsv_img(:,:,2) <= channel2Max) & ...
        (hsv_img(:,:,3) >= channel3Min ) & (hsv_img(:,:,3) <= channel3Max);
    mask = sliderBW;
    mask = imdilate(mask, strel('rectangle',[11 7])); 

   
end

