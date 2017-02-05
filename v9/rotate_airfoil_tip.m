function [X, R, Y] = rotate_airfoil_tip(X, R, Y, skew, reversex, phiplus)
    
    phi = (abs(atand(skew/R(1)))+phiplus)*reversex;
    
    RM = [cosd(phi) -sind(phi) 0;
          sind(phi) cosd(phi) 0;
          0 0 1];
    Rst = R;
    R = 0*R;
    V = RM*[X; R; Y];
    X = V(1,:);
    R = V(2,:);
    Y = V(3,:);
    R = R + Rst;
end