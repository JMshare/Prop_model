clear all;
clc;
close all;
format;


%% Inputs
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
Nx = 50;
if Nx<length(A(1))
    sprintf('Warning: Nx less than airfoil data.\n');
end
if mod(Nx,2)~=0
    sprintf('Warning: inserted odd Nx. Corrected to Nx-1.\n');
    Nx = Nx-1;
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



%% Transform, store, export and plot the airfoils
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




%%% Meshing
F = zeros(0,3);
V = zeros(0,3);


%% Mesh Blade 
plti = 1;
closetrail = 1;
[F, V, eg_st, eg_nd] = mesh_surface(Xstor, Rstor, Ystor, Nr, F, V, closetrail, plti, 8, 11);


% close the mesh on the tip 
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




%% Mesh Hub 
nBlades = 2;
alpha0 = 0;
totalpha = 360/nBlades;
alphaT = totalpha+alpha0;
Nu = 25;
Nd = 15;
Nb = 15;
Nf = 15;
Nh = 15;
[Fh, Vh] = hub_mesh(Xstor, Rstor, Ystor, HubD, HubT, ShaftD, alpha0, alphaT, Nu, Nd, Nb, Nf, Nh);




%% Connect Hub and Blade
F = [F; Fh+length(V)];
V = [V; Vh];




%% Postproc and outputs

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

fprintf('Generated %d Verticies and %d Faces.\n', length(V), length(F));

    
