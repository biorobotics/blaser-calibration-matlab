function radialDistortion (k1, k2, k3)
    max = 10;
    [x, y] = meshgrid(-max:.5:max);
    r2 = x.^2 + y.^2;
    k = k1*r2 + k2*r2.^2 + k3*r2.^3;
    figure;
    quiver(x, y, x.*k, y.*k, 0);
    axis("square");
    figure; 
    
end




