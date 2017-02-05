function [F, V] = hub_mesh(Xstor, Rstor, Ystor, HubD, HubT, ShaftD, alpha0, alphaT, Nu, Nd, Nb, Nf, Nh)
    

    F = zeros(0,3);
    V = zeros(0,3);

    % up from root
    [~, ixmn] = min(Xstor(1,:));
    eg_from = [Xstor(1,1:ixmn); Rstor(1,1:ixmn); Ystor(1,1:ixmn)];
    eg_to = [Xstor(1,1:ixmn); Rstor(1,1:ixmn); Ystor(1,1:ixmn)*0+HubT/2];
    eg_st = zeros(3,0); % store the start mesh segment edge
    eg_nd = zeros(3,0); % store the end mesh segment edge
    N = length(eg_from);
    [F, V, eg_st, eg_nd] = mesh_area(eg_from, eg_to, Nu, F, V);
    STup = eg_to;
    S1up = eg_st;
    S2up = eg_nd;

    % down from root
    eg_to = fliplr([Xstor(1,ixmn:end); Rstor(1,ixmn:end); Ystor(1,ixmn:end)]);
    eg_from = fliplr([Xstor(1,ixmn:end); Rstor(1,ixmn:end); Ystor(1,ixmn:end)*0-HubT/2]);
    [F, V, eg_st, eg_nd] = mesh_area(eg_from, eg_to, Nd, F, V);
    STdw = eg_from;
    S1dw = eg_st;
    S2dw = eg_nd;


    % back from root
    eg_from = [S2dw, S2up(:,2:end)];
    N = length(eg_from)
    ys = linspace(-HubT/2, HubT/2, N);
    rs = 0*ys + HubD*sind(alphaT)/2;
    xs = 0*ys + HubD*cosd(alphaT)/2;
    eg_to = [xs; rs; ys];
    [F, V, eg_st, eg_nd] = mesh_area(eg_from, eg_to, Nb, F, V);
    STdw = [STdw, eg_st(:,2:end)];
    STup = [STup, eg_nd(:,2:end)];


    % front from root
    eg_to = [S1dw, S1up];
    N = length(eg_to)
    ys = linspace(-HubT/2, HubT/2, N);
    rs = 0*ys + HubD*sind(alpha0)/2;
    xs = 0*ys + HubD*cosd(alpha0)/2;
    eg_from = [xs; rs; ys];
    [F, V, eg_st, eg_nd] = mesh_area(eg_from, eg_to, Nf, F, V);
    STdw = [eg_st, STdw(:,2:end)];
    STup = [eg_nd, STup(:,2:end)];


    % upper hub
    eg_from = STup;
    N = length(eg_from);
    phi = linspace(alpha0, alphaT, N);
    eg_to = [ShaftD*cosd(phi)/2; ShaftD*sind(phi)/2; STup(3,:)]; 
    [F, V, eg_st, eg_nd] = mesh_area(eg_from, eg_to, Nh, F, V);
    SSup = eg_to;

    % down hub
    eg_from = STdw;
    N = length(eg_from);
    phi = linspace(alpha0, alphaT, N);
    eg_to = [ShaftD*cosd(phi)/2; ShaftD*sind(phi)/2; STdw(3,:)]; 
    [F, V, eg_st, eg_nd] = mesh_area(eg_from, eg_to, Nh, F, V);
    SSdw = eg_to;

    % inner shaft
    eg_from = SSdw;
    eg_to = SSup;
    n_from = length(eg_from)
    n_to = length(eg_to)
    N = max(n_from, n_to);
    n_diff = abs(n_from-n_to);
    if n_diff~=0
        for i = 1:N-1
            if mod(i,floor((N-2)/n_diff)) == 0 % to fill in the missing points
                if n_from > n_to
                    eg_to = [eg_to(:,1:i), eg_to(:,i:end)];
                elseif n_from < n_to
                    eg_from = [eg_from(:,1:i), eg_from(:,i:end)]; % (the ith point repeats here)
                end
            end
        end
    end
    [F, V, eg_st, eg_nd] = mesh_area(eg_from, eg_to, Nu+Nd, F, V);


end