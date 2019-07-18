function tangentialDistortion (p1, p2)
    max = 10;
    [x, y] = meshgrid(-max:.5:max);
    r = x.^2 + y.^2;
    xy = x.*y*2;
    figure;
    quiver(x, y, xy*p1 + (r+2*x.^2)*p2, xy*p2  + (r+2*y.^2)*p1, 0);
    axis("square");
end