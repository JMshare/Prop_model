function [X, R, Y] = rotate_blade(X, R, Y, phi)


    RM = [cosd(phi) -sind(phi) 0;
          sind(phi) cosd(phi) 0;
          0 0 1];
    V = RM*[X; R; Y];
    X = V(1,:);
    R = V(2,:);
    Y = V(3,:);
end