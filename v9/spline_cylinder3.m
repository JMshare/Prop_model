function [S] = spline_cylinder3(points, N)
    
    RF = sqrt(points(1,1)^2+points(1,2)^2);
    RT = sqrt(points(2,1)^2+points(2,2)^2);
    Rs = linspace(RF, RT, N);
    
    phiF = acosd(points(1,1)/RF);
    phiT = acosd(points(2,1)/RT);
    phi = linspace(phiF, phiT, N);
    
    
    s1x = Rs.*cosd(phi);
    s1r = Rs.*sind(phi);
    s1y = linspace(points(1,3), points(2,3), N);
    
    S = [s1x; s1r; s1y];

end