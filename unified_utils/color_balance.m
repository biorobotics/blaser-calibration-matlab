function [out_im] = color_balance(img, percent, debug)
    img = double(img)/255;
    satLevel = percent/100;
    q = [satLevel/2 1-satLevel/2];
    imRGB_orig = cbreshape(img)*255;
    imRGB = zeros(size(imRGB_orig));
    N = size(imRGB_orig,2);
    color = {'r','g','b'};
    for ch = 1:3
        if debug
            figure;
            subplot(211);
            histogram(imRGB_orig(ch,:),256);
            set(findobj(gca,'Type','patch'),'FaceColor',color{ch},'EdgeColor',color{ch});
            xlim([0 255]);
            title('Original Histogram');
        end
        tiles = quantile(imRGB_orig(ch,:),q);
    %     [sum(imRGB_orig(ch,:)<tiles(1))/N,sum(imRGB_orig(ch,:)>tiles(2))/N] %check percentages are correct
        imRGB(ch,:) = cbsaturate(imRGB_orig(ch,:),tiles); %saturate at the appropriate pts. in distribution
        bottom = min(imRGB(ch,:)); top = max(imRGB(ch,:));
        imRGB(ch,:) = (imRGB(ch,:)-bottom)*255/(top-bottom);

        if debug
            subplot(212);
            histogram(imRGB(ch,:),256);
            set(findobj(gca,'Type','patch'),'FaceColor',color{ch},'EdgeColor',color{ch});
            xlim([0 255]);
            title('Corrected Histogram');
        end
    end


    out_im =cbunshape(imRGB,size(img))/255;
%     figure;
%     imshow(out_im);
%     title('Simplest Color Balance Corrected');
end

function im = cbimread(filename)
%im = cbimread(filename)
% Reads an image and returns it in 0-1 float format.

im = double(imread(filename))/255;
end

function out = cbreshape(im)
%out = cbreshape(im)
% Takes a width x height x 3 RGB image and returns a matrix where each column is an RGB
% pixel.

s = size(im);
m = s(1); n = s(2);
out = reshape(permute(im,[3 1 2]),[3 m*n 1]);

% x=[]
% for i=1:3
%     x=[x;reshape(im(:,:,i),m*n,1)'];
% end
end

function y = cbsaturate(x,bounds)
%y = cbsaturate(x,bounds)
% Saturates the vector x given the bounds pair [low high]

low = bounds(1); high = bounds(2);

y=zeros(size(x));
for i = 1:length(x)
    if x(i) < low
        y(i) = low;
    elseif x(i) > high
        y(i) = high;
    else
        y(i) = x(i);
    end
end
end

function out = cbunshape(mat,s)
% Takes a 3xn matrix of RGB pixels and returns a height x width x 3 RGB
% image

width = s(1); height = s(2);
out = reshape(mat,[3,width,height]);
out = permute(out,[2 3 1]);
end