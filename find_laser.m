function [pts] = find_laser(img, thresh)
    red = img(:,:,1);    
    pts = [];
    
    % Slice img along colums
    for i=1:size(img, 2)
        col = red(:,i);
        
        [vmax, imax] = max(col);
        
        if vmax <= thresh
            continue
        end
        
        [lmin, li] = min(col(1:imax));
        [rmin, ri] = min(col(imax:size(col, 1)));
        
        prom = vmax - lmin;
        if rmin > lmin
            prom = vmax - rmin;
        end

        vt = cast(vmax - 0.5*prom, 'double');
        xs = find(col > vt);
        
        % find discontinuous indices
        seps = [1, find(diff(xs) > 5)'+1, size(xs,1)+1];
        sex = find(seps > find(xs==imax), 1)-1;
        xs = xs(seps(sex) : seps(sex+1)-1);
        
        if size(xs,1) < 3
            continue
        end
        roi = cast(col(xs), 'double');
        % center of mass method
        ms = roi - vt;
        peak_loc = (ms'*xs)/(sum(ms));
        pts = [pts; i, peak_loc];
    end
end

% function [pts, pts_max] = find_laser(img)
%     red = img(:,:,1);
%     RED_THRESH = 200;
%     peakfind_percentile = 0.5;
%     
%     pts = [];
%     pts_max = [];
%     
%     % Slice img along colums
%     figure;
%     for i=1:size(img, 2)
%         col = red(:,i);
%         clf;
%         plot(col);
%         ylim([0,255]);
%         hold on;
%         
%         [pks,locs, ~,proms] = findpeaks(cast(col, 'single'),'MinPeakProminence',50);
%         [~,poi] = max(pks);
%         if size(pks, 1) == 0 || pks(poi) < RED_THRESH
%             drawnow;
%             continue
%         end
%         pts_max = [pts_max; i,locs(poi)];
%         scatter(locs(poi), pks(poi), 'go');
% %         drawnow;
% 
%         vt = pks(poi) - (1-peakfind_percentile)*proms(poi);
%         xs = find(col > vt);
%         
%         disp(i);
%         % find discontinuous indices
%         seps = [1, find(diff(xs) ~= 1)'+1, size(xs,1)+1];
%         sex = find(seps >= find(xs==locs(poi)), 1)-1;
%         xs = xs(seps(sex) : seps(sex+1)-1);
%         
%         if size(xs,1) < 3
%             disp('Too msall!"');
%             drawnow;
%             continue
%         end
%         roi = cast(col(xs), 'double');
%         % center of mass method
%         ms = roi - vt;
%         peak_loc = (ms'*xs)/(sum(ms));
%         line([peak_loc peak_loc],get(gca,'YLim'),'Color','r');
% %         vline(peak_loc,'r');
% 
%         % 2nd order poly fit
% %         obj = fit(xs, roi, 'poly2');
% %         peak_loc = -0.5*obj.p2/obj.p1;
% %         scatter(peak_loc, obj.p1*peak_loc*peak_loc + obj.p2*peak_loc + obj.p3, 'r+');
%         drawnow;
% 
% %         pause(0.3);
%         pts = [pts; i, peak_loc];
%     end
% end
