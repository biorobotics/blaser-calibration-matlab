function [] = show_reproj_err(A, ext, state)
% Shows reprojection errors
cm = hsv(54);

err_means = zeros(54,2);
nim = 0;

figure;
for i=1:size(A,1)
    if size(A{i,2}, 1) == 0
        continue
    end
    
    err = reproj_err(A{i,2}, ext, state);
    scatter(err(1:54), err(55:108),'Marker', '.','CData',cm);
    hold on;
    
    err_means(:,1) = err_means(:,1) + err(1:54);
    err_means(:,2) = err_means(:,2) + err(55:108);
    nim = nim+1;
end
scatter(err_means(:,1)/nim, err_means(:,2)/nim,'CData',cm, 'Marker','o');

end
