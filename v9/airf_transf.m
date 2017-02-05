function [X, Y] = airf_transf(A, cl, c, beta, rake, skew, thickness, reversex, reversey, Nx)

    % airfoil coordinates
    X = A(:,1);
    Y = A(:,end);
    if((X(end) == X(1)) && (Y(end) == Y(1)))
        X = X(1:end-1);
        Y = Y(1:end-1);
    end
    
    % chordline coordinates
    Xc = linspace(min(X),max(X), 100)';
    Yc = 0*Xc;
    
    % plot original
    figure(1)
    plot(X, Y);
    axis equal;
    grid on;
    set(gcf, 'color', 'white');
    title('Airfoil data original');
    hold on;
    plot(Xc, Yc, 'k--');
    
    % resize the airfoil, shift to the centerline and flip if needed.
    X = ((X - cl/100)*c)*reversex;
    Y = (Y*thickness*reversey/100)*c;
    Xc = ((Xc - cl/100)*c)*reversex;
    Yc = (Yc*thickness*reversey/100)*c;
 
    
    % rotate it
    beta = -1*beta*reversex;
    RM = [cosd(beta) -sind(beta);
          sind(beta) cosd(beta)];
    V = RM*[X,Y]';
    X = V(1,:);
    Y = V(2,:);
    Vc = RM*[Xc,Yc]';
    Xc = Vc(1,:);
    Yc = Vc(2,:);
  
    
    % shift it
    X = X+skew;
    Y = Y+rake;
    Xc = Xc+skew;
    Yc = Yc+rake;
    
    % plot transformed
    figure(2)
    plot(X, Y, 'g');
    hold on;
    plot(Xc, Yc, 'k--');
    title('Transformed airfoil');
    grid on;
    
    
    % Spline interpolate it
    lng = length(X);
    points = [X; Y];
    pp = cscvn(points);
    t = cumsum([0;((diff(points.').^2)*ones(2,1)).^(1/4)]).';
    splup = fnval(pp, linspace(t(1), t(floor(lng/2)), Nx/2));
    spldw = fnval(pp, linspace(t(floor(lng/2)), t(end), Nx/2+1));
    spl = [splup, spldw(:,2:end)];
    X = spl(1,:);
    Y = spl(2,:);
    
    
    figure(5)
    plot(X, Y, 'b.');
    axis equal;
    grid on;
    set(gcf, 'color', 'white');
    title('Airfoil splined');
    
end


