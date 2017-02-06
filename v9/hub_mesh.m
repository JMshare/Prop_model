function [F, V] = hub_mesh(Xstor, Rstor, Ystor, HubD, HubT, ShaftD, alpha0, alphaT, Nu, Nd, Nb, Nf, Nh)
    

    F = zeros(0,3);
    V = zeros(0,3);

    % up from root
    [~, ixmn] = min(Xstor(1,:));
    eg_from = [Xstor(1,1:ixmn); Rstor(1,1:ixmn); Ystor(1,1:ixmn)];
    eg_to = [Xstor(1,1:ixmn); Rstor(1,1:ixmn); Ystor(1,1:ixmn)*0+HubT/2];
    [F, V, eg_st, eg_nd] = mesh_area(eg_from, eg_to, Nu, F, V);
    STup = eg_to;
    S1up = eg_st;
    S2up = eg_nd;

    % down from root
    eg_to = fliplr([Xstor(1,ixmn:end); Rstor(1,ixmn:end); Ystor(1,ixmn:end)]);
    eg_from = fliplr([Xstor(1,ixmn+1:end); Rstor(1,ixmn+1:end); Ystor(1,ixmn+1:end)*0-HubT/2]);
    [F, V, eg_st, eg_nd] = mesh_area(eg_from, eg_to, Nd, F, V);
    STdw = eg_from;
    S1dw = eg_st;
    S2dw = eg_nd;


    % back from root
    eg_from = [S2dw, S2up];
    ys = linspace(-HubT/2, HubT/2, Nu+Nd);
    rs = 0*ys + HubD*sind(alphaT)/2;
    xs = 0*ys + HubD*cosd(alphaT)/2;
    eg_to = [xs; rs; ys];
    [eg_from, eg_to, N] = blend_edges(eg_from, eg_to, 3);
    N
    [F, V, eg_st, eg_nd] = mesh_area(eg_from, eg_to, Nb, F, V);
    STdw = [STdw, eg_st(:,2:end)];
    STup = [STup, eg_nd(:,2:end)];


    % front from root
    eg_to = [S1dw, S1up];
    ys = linspace(-HubT/2, HubT/2, Nu+Nd);
    rs = 0*ys + HubD*sind(alpha0)/2;
    xs = 0*ys + HubD*cosd(alpha0)/2;
    eg_from = [xs; rs; ys];
    [eg_from, eg_to, N] = blend_edges(eg_from, eg_to, 4);
    N
    [F, V, eg_st, eg_nd] = mesh_area(eg_from, eg_to, Nf, F, V);
    STdw = [eg_st, STdw(:,2:end)];
    STup = [eg_nd, STup(:,2:end)];


    % upper hub
    eg_from = STup;
    N = length(eg_from);
    phi = linspace(alpha0, alphaT, N);
    eg_to = [ShaftD*cosd(phi)/2; ShaftD*sind(phi)/2; STup(3,:)]; 
    [F, V, ~, ~] = mesh_area(eg_from, eg_to, Nh, F, V);
    SSup = eg_to;

    % down hub
    eg_from = STdw;
    N = length(eg_from);
    phi = linspace(alpha0, alphaT, N);
    eg_to = [ShaftD*cosd(phi)/2; ShaftD*sind(phi)/2; STdw(3,:)];
    [F, V, ~, ~] = mesh_area(eg_from, eg_to, Nh, F, V);
    SSdw = eg_to;

    % inner shaft
    eg_from = SSdw;
    eg_to = SSup;
    [eg_from, eg_to, ~] = blend_edges(eg_from, eg_to, 7);
    [F, V, ~, ~] = mesh_area(eg_from, eg_to, Nu+Nd, F, V);


end