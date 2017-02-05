function [S] = spline_cubic3(points, N)

    % generate cubic spline 
    pp = cscvn(points);
    t = cumsum([0;((diff(points.').^2)*ones(3,1)).^(1/4)]).';
    S = fnval(pp, linspace(min(t), max(t), N));

end