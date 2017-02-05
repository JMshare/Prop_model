clear all;
clc;
close all;
format;

R = 9*25.4/2; % [mm] Radius of the propeller
HubD = 0.5*25.4;
HubT = 0.31*25.4;
ShaftD = 0.25*25.4;
cl = 30; % [%] position of the centerline (origin) to the chord (50% is mid chord)
thicknessm = 100; % [%] thickness of the airfoil compared to the original one
reversey = 1; % [1/-1] set -1 to reverse the airfoil vertically
reversex = 1; % [1/-1] set -1 to reverse the airfoil horizontally

C = dlmread('apcsf_9x4.7_geom.txt', ' ', 1, 0); % from http://m-selig.ae.illinois.edu/props/volume-1/propDB-volume-1.html. Make sure the format is all right.
bin = not(all(C==0,1));
C = C(:,bin);
r = C(:,1)*R;
c = C(:,2)*R;
beta = C(:,3);
rake = -r*HubT/r(end) + HubT; % vertical shift
skewm = HubD/2; % horizontal shift
skewshft = 3;
skew = sin(r*pi/R)*skewm + skewshft;
thickness = r*0.0 + thicknessm;

c(1) = c(1)/1.;



A = {};
for i=1:length(r)
    D = importdata('naca4412.dat', ' ', 1); % from http://m-selig.ae.illinois.edu/ads/coord_database.html. Make sure the format is all right. Must start from X=1, not 0
    %D = importdata('clarkysm.dat');
    A(i) = {D.data};
end
A = A';


Nr = 100;
Nx = 26;
if Nx<length(A(1))
    sprintf('Warning: Nx less than airfoil data.\n');
end


% Tip airfoils
phiTip = 20;
dTip = 0.5;
rTip = r(end)+dTip;
cTip = c(end)/2.0;
if c(end)/cTip < 0.2
    sprintf('Warning, tip too squeezed. Consider rotating it.\n');
end
Y = cell2mat(A(end));
Y = Y(:,2);
tTip = min(abs(Y(abs(Y)>0)));
tTip = tTip + min(abs(Y(abs(Y)>tTip)));
skewTip = -1.8+skewshft;
Nt = 8;
[r, c, beta, rake, skew, thickness, A] = generate_Tip(r, c, beta, rake, skew, thickness, A, dTip, rTip, cTip, tTip, skewTip, Nt);


% Root airfoils
dRoot = r(1)-HubD/2;
rRoot = HubD/2;
cRoot = pi*HubD/3.2;
if cRoot > HubD
    sprintf('Warning, root chord cant exceed the hub diameter.\n');
end
tRoot = (HubT/3)/HubD;
betaRoot = 0;
skewRoot = 2;
rakeRoot = 1;
thickRoot = 130;
Nrt = 3;
[r, c, beta, rake, skew, thickness, A] = generate_Root(r, c, beta, rake, skew, thickness, A, dRoot, rRoot, cRoot, tRoot, betaRoot, skewRoot, thickRoot, rakeRoot, Nrt);



% Transform, store, export and plot the airfoils
figure(3)
hold on;
Xstor = zeros(length(r), Nx);
Ystor = Xstor;
Rstor = Xstor;
for i = 1:length(r) % for each airfoil along r
    % transform airfoil
    [X, Y] = airf_transf(cell2mat(A(i)), cl, c(i), beta(i), rake(i), 0*skew(i), thickness(i), reversex, reversey, Nx);
    R = r(i)+X*0;
    
    % project airfoil to a cylinder
    [X, R, Y] = project_airfoil(X, R, Y, r(i));
    
    % rotate the airfoils around the center line to create the skew
    [X, R, Y] = rotate_airfoil(X, R, Y, r(i), skew(i), reversex, 0);
    
    % rotate tip blades
    if i >= (length(r) - Nt)
        [X, R, Y] = rotate_airfoil_tip(X, R, Y, 0, reversex, phiTip);
    end 
    
    % store airfoil
    Xstor(i,:) = X;
    Ystor(i,:) = Y;
    Rstor(i,:) = R;
   
%     % Write to solidworks text file
%     str = sprintf('./Matlab_Out/af%d.sldcrv', i);
%     fileID = fopen(str, 'w');
%     nbytes = fprintf(fileID,'%f %f %f\r\n', [[X]; [R]; [Y]]);
%     fclose(fileID);
% 
%     % Write to .dxf
%     str = sprintf('./Matlab_Out/af%d.dxf', i);
%     FID = dxf_open(str);
%     dxf_polyline(FID, [X, X(1)], [R, R(1)], [Y, Y(1)]);
%     dxf_close(FID);
%     
%     % Write to iges
%     str = sprintf('./Matlab_Out/af%d.iges', i);
%     igesout([[X, X(1)]', [R, R(1)]', [Y, Y(1)]'], str);

end

% plot the skeleton of transformed airfoils and hub
plot_skeleton(3, Xstor, Rstor, Ystor, r, skew, rake, HubD, HubT, ShaftD);

% plot tip
figure(6)
hold on;
for i = (length(r)-Nt): length(r)
    plot3(Xstor(i,:),  Rstor(i,:), Ystor(i,:), 'b');
end
axis equal;
title('Tip profile');

% plot tip airfoil separately
figure(9)
plot(Xstor(end,:), Ystor(end,:), 'b');
title('Tip airfoil separately');

% Plot cubic splines along r
figure(4)
hold on;
for i = 1:Nx
    fnplt(cscvn([Xstor(:,i)'; Rstor(:,i)'; Ystor(:,i)']),'r',2);
end
axis equal;



%%% Meshing
F = zeros(0,3);
V = zeros(0,3);


% Create mesh along r 
plti = 1;
closed = 1;
N = size(Xstor, 2); % Nx
for i=1:N-(closed==0)
    [S1] = spline_cubic3([Xstor(:,i)'; Rstor(:,i)'; Ystor(:,i)'], Nr);
    i = mod(i,N); % to get 0 if i==Nx so that it can close
    [S2] = spline_cubic3([Xstor(:,i+1)'; Rstor(:,i+1)'; Ystor(:,i+1)'], Nr);
    
    [Fs, Vs, F, V] = generate_mesh_slice(F, V, S1, S2);
    
    Nr = length(S2);
    
    if i == plti
        figure(8)
        axis equal;
        hold on;
        trimesh(Fs, Vs(:,1),Vs(:,2), Vs(:,3));
        title('The whole observed slice');
        figure(11)
        trimesh(Fs(end-3:end,:), Vs(:,1),Vs(:,2), Vs(:,3), 'FaceColor', 'blue');
        axis equal;
        hold on;
    end
    
end

% close the mesh on the tip mesh
X = Xstor(end,:);
Y = Ystor(end,:);
R = Rstor(end,:);
S2 = [X(1:Nx/2); R(1:Nx/2); Y(1:Nx/2)];
S1 = fliplr([X(Nx/2+1:end); R(Nx/2+1:end); Y(Nx/2+1:end)]);
[Fs, Vs, F, V] = generate_mesh_slice(F, V, S1, S2);
figure(10)
trimesh(Fs, Vs(:,1),Vs(:,2), Vs(:,3));
axis equal;
title('The tip edge');
view(0,45);
figure(11);
hold on;
trimesh(Fs(1:3, :), Vs(:,1),Vs(:,2), Vs(:,3), 'FaceColor', 'red');
axis equal;
title('The tip edge bit with the trailing edge bit');





% create Hub mesh
nBlades = 2;
totalpha = 360/nBlades;
alpha0 = 0;

Nh = 25;
[~, ixmn] = min(Xstor(1,:));
eg_from = [Xstor(1,1:ixmn); Rstor(1,1:ixmn); Ystor(1,1:ixmn)];
eg_to = [Xstor(1,1:ixmn); Rstor(1,1:ixmn); Ystor(1,1:ixmn)*0+HubT/2];
eg_st = zeros(3,0); % store the start mesh segment edge
eg_nd = zeros(3,0); % store the end mesh segment edge
N = length(eg_from);
for i = 1:N-1 % up from root
    S1 = spline_cylinder3([eg_from(:,i)'; eg_to(:,i)'],...
                     Nh);
    S2 = spline_cylinder3([eg_from(:,i+1)'; eg_to(:,i+1)'],... 
                     Nh - 1*(eg_from(3,i+1)>eg_from(3,i)) + 1*(eg_from(3,i+1)<eg_from(3,i)));
    
    [~, ~, F, V] = generate_mesh_slice(F, V, S1, S2);
    
    Nh = length(S2);
    
    
    if i == 1
        eg_st = S1;
    end
    if i == N-1
        eg_nd = S2;
    end
            
end
STup = eg_to;
S1up = eg_st;
S2up = eg_nd;


Nh = 15;
eg_to = fliplr([Xstor(1,ixmn:end); Rstor(1,ixmn:end); Ystor(1,ixmn:end)]);
eg_from = fliplr([Xstor(1,ixmn:end); Rstor(1,ixmn:end); Ystor(1,ixmn:end)*0-HubT/2]);
eg_st = zeros(3,0); % store the start mesh segment edge
eg_nd = zeros(3,0); % store the end mesh segment edge
N = length(eg_from);
for i = 1:N-1 % down from rooot
    S1 = spline_cylinder3([eg_from(:,i)'; eg_to(:,i)'],...
                     Nh);
    S2 = spline_cylinder3([eg_from(:,i+1)'; eg_to(:,i+1)'],... 
                     Nh + 1*(eg_to(3,i+1)>eg_to(3,i)) - 1*(eg_to(3,i+1)<eg_to(3,i)));
    
    [~, ~, F, V] = generate_mesh_slice(F, V, S1, S2);
    
    Nh = length(S2);
    
    if i == 1
        eg_st = S1;
    end
    if i == N-1
        eg_nd = S2;
    end
end
STdw = eg_from;
S1dw = eg_st;
S2dw = eg_nd;



alphaT = totalpha+alpha0;
Nb = 15;
eg_from = [S2dw, S2up(:,2:end)];
N = length(eg_from)
ys = linspace(-HubT/2, HubT/2, N);
rs = 0*ys + HubD*sind(alphaT)/2;
xs = 0*ys + HubD*cosd(alphaT)/2;
eg_to = [xs; rs; ys];
eg_st = zeros(3,0); % store the start mesh segment edge
eg_nd = zeros(3,0); % store the end mesh segment edge
for i = 1:N-1 % back from root
    S1 = spline_cylinder3([eg_from(:,i)'; eg_to(:,i)'], Nb);
    S2 = spline_cylinder3([eg_from(:,i+1)'; eg_to(:,i+1)'], Nb);
    
    [~, ~, F, V] = generate_mesh_slice(F, V, S1, S2);
    
    Nb = length(S2);
    
    if i == 1
        eg_st = S1;
    elseif i == N-1
        eg_nd = S2;
    end
end
STdw = [STdw, eg_st(:,2:end)];
STup = [STup, eg_nd(:,2:end)];



Nf = 15;
eg_to = [S1dw, S1up];
N = length(eg_to)
ys = linspace(-HubT/2, HubT/2, N);
rs = 0*ys + HubD*sind(alpha0)/2;
xs = 0*ys + HubD*cosd(alpha0)/2;
eg_from = [xs; rs; ys];
for i = 1: N-1 % front from root
    S1 = spline_cylinder3([eg_from(:,i)'; eg_to(:,i)'], Nf);
    S2 = spline_cylinder3([eg_from(:,i+1)'; eg_to(:,i+1)'], Nf);
    
    [~, ~, F, V] = generate_mesh_slice(F, V, S1, S2);
    
    Nf = length(S2);
    
    if i == 1
        eg_st = S1;
    elseif i == N-1
        eg_nd = S2;
    end
    
end
STdw = [eg_st, STdw(:,2:end)];
STup = [eg_nd, STup(:,2:end)];


figure(100)
plot3(STup(1,:), STup(2,:), STup(3,:));
hold on;
plot3(STdw(1,:), STdw(2,:), STdw(3,:));
axis equal;


Nu = 15;
eg_from = STup;
N = length(eg_from);
phi = linspace(alpha0, alphaT, N);
eg_to = [ShaftD*cosd(phi)/2; ShaftD*sind(phi)/2; STup(3,:)]; 
for i = 1: N-1 % upper Hub
    S1 = spline_cylinder3([eg_from(:,i)'; eg_to(:,i)'], Nu);
    S2 = spline_cylinder3([eg_from(:,i+1)'; eg_to(:,i+1)'], Nu);
    
    [~, ~, F, V] = generate_mesh_slice(F, V, S1, S2);
    
    Nu = length(S2);
    
    if i == 1
        eg_st = S1;
    elseif i == N-1
        eg_nd = S2;
    end     
    
end
SSup = eg_to;


Nu = 15;
eg_from = STdw;
N = length(eg_from);
phi = linspace(alpha0, alphaT, N);
eg_to = [ShaftD*cosd(phi)/2; ShaftD*sind(phi)/2; STdw(3,:)]; 
for i = 1: N-1 % upper Hub
    S1 = spline_cylinder3([eg_from(:,i)'; eg_to(:,i)'], Nu);
    S2 = spline_cylinder3([eg_from(:,i+1)'; eg_to(:,i+1)'], Nu);
    
    [~, ~, F, V] = generate_mesh_slice(F, V, S1, S2);
    
    Nu = length(S2);
    
    if i == 1
        eg_st = S1;
    elseif i == N-1
        eg_nd = S2;
    end
    
end
SSdw = eg_to;


Nh = 25;
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

for i = 1: N-1 % upper Hub
    S1 = spline_cylinder3([eg_from(:,i)'; eg_to(:,i)'], Nh);
    S2 = spline_cylinder3([eg_from(:,i+1)'; eg_to(:,i+1)'], Nh);
    
    [~, ~, F, V] = generate_mesh_slice(F, V, S1, S2);
    
    Nh = length(S2);
    
    if i == 1
        eg_st = S1;
    elseif i == N-1
        eg_nd = S2;
    end
    
end



figure(7)
plot_skeleton(7, Xstor, Rstor, Ystor, r, skew, rake, HubD, HubT, ShaftD);
trimesh(F, V(:,1),V(:,2), V(:,3));
axis equal;


% Rotate the whole part to get N blades
[X2, R2, Y2] = rotate_blade(V(:,1)',V(:,2)', V(:,3)', totalpha);
Vs = [X2', R2', Y2'];
Fs = F;
F = [F; Fs+length(V)];
V = [V; Vs];
figure(13)
trimesh(F, V(:,1),V(:,2), V(:,3));
axis equal;
plot_skeleton(13, Xstor, Rstor, Ystor, r, skew, rake, HubD, HubT, ShaftD);



% Remove duplicate vertices (the stl file will add them back again anyway though)
[V, indexm, indexn] =  unique(V, 'rows', 'stable');
F = indexn(F);

% Write to .stl CAD model
stlwrite('./Matlab_Out/test.stl', F, V, 'mode', 'ascii');

% Write to .dxf
% FID = dxf_open('./Matlab_Out/test1.dxf');
% dxf_polymesh(FID, V, F);
% dxf_close(FID);


figure(7)
    
