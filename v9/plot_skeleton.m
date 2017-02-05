function [] = plot_skeleton(fign, Xstor, Rstor, Ystor, r, skew, rake, HubD, HubT, ShaftD)

    figure(fign)
    hold on;

    % plot airfoils
    for i = 1:length(r)
        plot3(Xstor(i,:), Rstor(i,:), Ystor(i,:), 'b');
    end
    axis equal;

    % plot centerline
    plot3(0*r, r, 0*r, 'k--');
    plot3(-skew, r, rake, 'b--');

    % plot Hub
    ang=0:0.01:2*pi; 
    xp=(HubD/2)*cos(ang);
    rp=(HubD/2)*sin(ang);
    yp = xp*0-HubT/2;
    plot3(xp, rp, yp, 'b');
    yp = xp*0+HubT/2;
    plot3(xp, rp, yp, 'b');
    xp=(ShaftD/2)*cos(ang);
    rp=(ShaftD/2)*sin(ang);
    yp = xp*0-HubT/2;
    plot3(xp, rp, yp, 'b');
    yp = xp*0+HubT/2;
    plot3(xp, rp, yp, 'b');

end