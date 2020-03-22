function [pts] = find_laser_new(img)

    pts = [];

    % Median Blur
    img = medfilt3(img);
    ksize = [7 7];
    sigma = [0.3*((ksize(2)-1)*0.5 - 1) + 0.8, ...
             0.3*((ksize(1)-1)*0.5 - 1) + 0.8];
    
    % Gaussian Blur
    img = imgaussfilt(img,sigma, 'FilterSize',ksize);
    
    % Color Balance
    %im_correct = color_balance(img,0.1, 0);
    
    % Red and Saturation masks
    rmask = getredmask(img);
    smask = getsatmask(img);
    
    %se = strel('rectangle',[30,10]);
    %extended_smask = imdilate(smask, se);
    
    %h = imshow(extended_smask);
    %set(h, 'AlphaData', 0.2);
    
    fprintf("smask=%.4f, rmask=%.4f\n", ...
        sum(smask, 'all'), sum(rmask, 'all'));
    
    % Different masking mode based on whether laser line is mostly
    % saturated.
    %if sum(extended_smask, 'all') > sum(rmask, 'all')
    if sum(smask, 'all') > size(smask, 2) * 2
        % laser stripe mostly saturated, search in saturation
        fprintf("use saturation mask instead\n");
        tmask = smask;
    else
        % laser stripe mostly not saturated, subtract saturation
        % when subtrarcting, make sure subtraction is thorough
        se = strel('rectangle',[300,30]);
        extended_smask = imdilate(smask, se);
        tmask = rmask & (~extended_smask);
    end
    im_comb = bitwise_and(img,tmask);
    
    % Fill 'air gaps' in mask
    se = strel('disk', 5);
    im_comb = imclose(im_comb, se);
    
    %h = imshow(rmask);
    %set(h, 'AlphaData', 0.5);
    
    %h = imshow(smask);
    %set(h, 'AlphaData', 0.7);
    
    m_img = mean(im_comb, 'all') / 2;
    im_sub_mean = im_comb - m_img;
    
    im_sub_mean(im_sub_mean < 0)=0;

    %imshow(im_comb)

    [val_max, col_max] = max(im_sub_mean, [], 1); 
    for i=1:size(img,2)
        j = col_max(i)-1;
   
        while (j >= 1 && im_sub_mean(j,i) > 0.8*val_max(i))
            j = j - 1;
        end
        k = col_max(i) + 1;
        
        while (k <= size(img,1) && im_sub_mean(k,i) > 0.8*val_max(i))
            k = k + 1;
        end
        weighted_sum = 0.0;
        val_sum = 0.0;
        for idx=j+1:k-1
            weighted_sum = weighted_sum + im_sub_mean(idx,i)*(idx-j);
            val_sum = val_sum + im_sub_mean(idx,i);
        end
        if (val_sum <=0)
            continue;
        end
        pts = [pts; i, double(weighted_sum/val_sum + j)];
    end
    
end


function [res] = bitwise_and(img,mask)

    res = zeros(size(img));
    for i = 1:size(img,1)
        for j = 1:size(img,2)
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
    channel1Min = 0.839;
    channel1Max = 0.050;
    % Define thresholds for channel 2 based on histogram settings
    channel2Min = 0.000;
    channel2Max = 1.000;
    % Define thresholds for channel 3 based on histogram settings
    channel3Min = 0.400;
    channel3Max = 1.000;
    % Create mask based on chosen histogram thresholds
    sliderBW = ( (hsv_img(:,:,1) >= channel1Min) | (hsv_img(:,:,1) <= channel1Max) ) & ...
        (hsv_img(:,:,2) >= channel2Min ) & (hsv_img(:,:,2) <= channel2Max) & ...
        (hsv_img(:,:,3) >= channel3Min ) & (hsv_img(:,:,3) <= channel3Max);
    mask = sliderBW;
    
end

function [mask] = getsatmask(img)

    % Convert RGB image to chosen color space
    hsv_img = rgb2hsv(img);
    % Define thresholds for channel 1 based on histogram settings
    channel1Min = 0.928;
    channel1Max = 0.100;
    % Define thresholds for channel 2 based on histogram settings
    channel2Min = 0.000;
    channel2Max = 0.209;
    % Define thresholds for channel 3 based on histogram settings
    channel3Min = 0.903;
    channel3Max = 1.000;
    % Create mask based on chosen histogram thresholds
    sliderBW = ( (hsv_img(:,:,1) >= channel1Min) | (hsv_img(:,:,1) <= channel1Max) ) & ...
        (hsv_img(:,:,2) >= channel2Min ) & (hsv_img(:,:,2) <= channel2Max) & ...
        (hsv_img(:,:,3) >= channel3Min ) & (hsv_img(:,:,3) <= channel3Max);
    mask = sliderBW;
    
    % Fill 'air gaps' in mask
    se = strel('disk', 5);
    mask = imclose(mask, se);
   
end

