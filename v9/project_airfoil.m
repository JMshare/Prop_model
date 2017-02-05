function [X, R, Y] = project_airfoil(X, R, Y, r)
    

%     alpha = abs(atan(R./X));
%     r = R(1);
%     X = r.*cos(alpha).*sign(X);
%     R = r.*sin(alpha);

    alpha = abs(X/r);
    R = r*cos(alpha);
    X = r*sin(alpha).*sign(X);
end