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


Nr = 200;
Nx = 50;
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




figure(3)
hold on;
Xstor = zeros(length(r-1), Nx);
Ystor = Xstor;
Rstor = Xstor;
for i = 1:length(r) % for each airfoil along r
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
    
    
    Xstor(i,:) = X;
    Ystor(i,:) = Y;
    Rstor(i,:) = R;
    
%     VV = fnplt(cscvn([Xstor(i,:); Rstor(i,:); Ystor(i,:)]));
%     str = sprintf('../v2/Matlab_Out/af%d.sldcrv', i);
%     fileID = fopen(str, 'w');
%     nbytes = fprintf(fileID,'%f %f %f\r\n', VV);
%     fclose(fileID);
end

plot_skeleton(3, Xstor, Rstor, Ystor, r, skew, rake, HubD, HubT, ShaftD);
plot_skeleton(7, Xstor, Rstor, Ystor, r, skew, rake, HubD, HubT, ShaftD);



% plot tip
figure(6)
hold on;
for i = (length(r)-Nt): length(r)
    plot3(Xstor(i,:),  Rstor(i,:), Ystor(i,:), 'b');
end
axis equal;
title('Tip');

% plot tip airfoil
figure(9)
plot(Xstor(end,:), Ystor(end,:), 'b');
title('Tip airfoil');


figure(4)
hold on;
% Create splines along r
Sstor = zeros(0, Nr);
for i = 1:Nx
    % plot just the sparse points
    fnplt(cscvn([Xstor(:,i)'; Rstor(:,i)'; Ystor(:,i)']),'r',2);
    
    % generate spline and fine points
    points = [Xstor(:,i)'; Rstor(:,i)'; Ystor(:,i)'];
    pp = cscvn(points);
    t = cumsum([0;((diff(points.').^2)*ones(3,1)).^(1/4)]).';
    spl = fnval(pp, linspace(min(t), max(t), Nr));
    
    % store the spline points for meshing later
    Sstor = [Sstor; spl];
end
axis equal;


% create the mesh slices and store them together
F = zeros(0,3);
V = zeros(0,3);
plti = 1;
for i=1:Nx
    ind1 = (i-1)*3+1 : i*3;
    if i == Nx
        i = 0;
    end
    ind2 = i*3+1 : (i+1)*3;
    S1 = Sstor(ind1, :);
    S2 = Sstor(ind2, :);
    % plot3(S1(1,:),S1(2,:),S1(3,:), 'b');
    % plot3(S2(1,:),S2(2,:),S2(3,:), 'b');
    S = [S1, S2];
    s1 = 1:Nr;
    s2 = Nr+1:2*Nr;
    s = [s1; s2];
    bin = s(:);
    Vs = S(:,bin)';
    Fs = [1:length(S)-2; 2:length(S)-1; 3:length(S)]';
    %patch('Faces', Fs, 'Vertices', Vs); % 'Faces' is 'Face' in
    % the documentation
    
    if i == plti
        figure(8)
        axis equal;
        hold on;
        trimesh(Fs, Vs(:,1),Vs(:,2), Vs(:,3));
        title('The whole plti slice');
        figure(11)
        trimesh(Fs(end-10:end,:), Vs(:,1),Vs(:,2), Vs(:,3), 'FaceColor', 'blue');
        axis equal;
        hold on;
    end
    
    F = [F; Fs+length(V)];
    V = [V; Vs];
end

% close the mesh on the tip
X = Xstor(end,:);
Y = Ystor(end,:);
R = Rstor(end,:);
S2 = [X(1:Nx/2); R(1:Nx/2); Y(1:Nx/2)];
S1 = fliplr([X(Nx/2+1:end); R(Nx/2+1:end); Y(Nx/2+1:end)]);
S = [S1, S2];
s1 = 1:length(S1);
s2 = length(S1)+1:length(S);
s = [s1; s2];
bin = s(:);
Vs = S(:,bin)';
Fs = [1:length(S)-2; 2:length(S)-1; 3:length(S)]';
figure(10)
hold on;
trimesh(Fs, Vs(:,1),Vs(:,2), Vs(:,3));
axis equal;
title('The tip edge');
view(0,45);
figure(11);
hold on;
trimesh(Fs(1:10, :), Vs(:,1),Vs(:,2), Vs(:,3), 'FaceColor', 'red');
axis equal;
title('The tip edge bit with the trailing edge bit');

F = [F; Fs+length(V)];
V = [V; Vs];


% create Hub
totalpha = 360/2;
alpha0 = 0;

[~, ixmn] = min(Xstor(1,:));
upairfx = Xstor(1,1:ixmn);
upairfy = Ystor(1,1:ixmn);
upairfr = Rstor(1,1:ixmn);
%dwairf = [Xstor(1,ixmn:end); Ystor(1,ixmn:end)];
Nh = 25;
STup = zeros(3,0); % store the upper hub edge
for i = 1:length(upairfx)-1 % up from root
    s1y = linspace(upairfy(i), HubT/2, Nh);
    s1x = s1y*0 + upairfx(i);
    s1r = s1y*0 + upairfr(i);
    s2y = linspace(upairfy(i+1), HubT/2, Nh - 1*(upairfy(i+1)>upairfy(i)) + 1*(upairfy(i+1)<upairfy(i)));
    s2x = s2y*0 + upairfx(i+1);
    s2r = s2y*0 + upairfr(i+1);
    S1 = [s1x; s1r; s1y];
    S2 = [s2x; s2r; s2y];
    if length(s1x) >= length(s2x)
        S = [S1, S2];
        S(:,1:2:length(S)) = S1;
        S(:,2:2:length(S)) = S2;
    else
        S = [S2, S1];
        S(:,1:2:length(S)) = S2;
        S(:,2:2:length(S)) = S1;
    end
    
    Vs = S';
    Fs = [1:length(S)-2; 2:length(S)-1; 3:length(S)]';
    Nh = length(s2x);
    
    F = [F; Fs+length(V)];
    V = [V; Vs];
    
    
    STup = [STup, S1(:,end)];
    if i == 1
        S1up = S1;
    end
    if i == length(upairfx)-1
        S2up = S2;
        STup = [STup, S2(:,end)];
    end
            
end

dwairf = [Xstor(1,ixmn:end); Rstor(1,ixmn:end); Ystor(1,ixmn:end)];
Nh = 15;
STdw = zeros(3,0); % store the down hub edge
for i = 1:length(dwairf)-1 % down from rooot
    s1y = linspace(dwairf(3,i), -HubT/2, Nh);
    s1x = s1y*0 + dwairf(1,i);
    s1r = s1y*0 + dwairf(2,i);
    s2y = linspace(dwairf(3,i+1), -HubT/2, Nh + 1*(dwairf(3,i+1)>dwairf(3,i)) - 1*(dwairf(3,i+1)<dwairf(3,i)));
    s2x = s2y*0 + dwairf(1,i+1);
    s2r = s2y*0 + dwairf(2,i+1);
    S1 = [s1x; s1r; s1y];
    S2 = [s2x; s2r; s2y];
    if length(s1x) >= length(s2x)
        S = [S1, S2];
        S(:,1:2:length(S)) = S1;
        S(:,2:2:length(S)) = S2;
    else
        S = [S2, S1];
        S(:,1:2:length(S)) = S2;
        S(:,2:2:length(S)) = S1;
    end
    
    Vs = S';
    Fs = [1:length(S)-2; 2:length(S)-1; 3:length(S)]';
    Nh = length(s2x);
    
    F = [F; Fs+length(V)];
    V = [V; Vs];
    
    STdw = [S1(:,end), STdw];
    if i == 1
        S1dw = S1;
    end
    if i == length(dwairf)-1
        S2dw = S2;
        STdw = [S2(:,end), STdw];
    end
end

alphaT = totalpha+alpha0;
Nb = 15;
edge1 = [fliplr(S1dw), S2up(:,2:end)];
ys = linspace(-HubT/2, HubT/2, length(edge1));
rs = 0*ys + HubD*sind(alphaT)/2;
xs = 0*ys + HubD*cosd(alphaT)/2;
edgeend = [xs; rs; ys];
phiF = acosd((edge1(1,1))/(HubD/2));
phi = linspace(phiF, alphaT, Nb);
length(edge1)
for i = 1:length(edge1)-1 % back from root
    
    xs = HubD*cosd(phi)/2;
    rs = HubD*sind(phi)/2;
    ys = linspace(edge1(3,i), edgeend(3,i), Nb);
    S1 = [xs; rs; ys];
    
    xs = HubD*cosd(phi)/2;
    rs = HubD*sind(phi)/2;
    ys = linspace(edge1(3,i+1), edgeend(3,i+1), Nb);
    S2 = [xs; rs; ys];
    
    S = [S1, S2];
    S(:,1:2:length(S)) = S1;
    S(:,2:2:length(S)) = S2;
    
    Vs = S';
    Fs = [1:length(S)-2; 2:length(S)-1; 3:length(S)]';
    
    F = [F; Fs+length(V)];
    V = [V; Vs];
    
    if i == 1
        STdw = [STdw, S1(:,2:end)];
    elseif i == length(edge1)-1
        STup = [STup, S2(:,2:end)];
    end
end



Nf = 15;
edge1 = [fliplr(S2dw), S1up];
ys = linspace(-HubT/2, HubT/2, length(edge1));
rs = 0*ys + HubD*sind(alpha0)/2;
xs = 0*ys + HubD*cosd(alpha0)/2;
edgeend = [xs; rs; ys];
phiF = acosd((edge1(1,1))/(HubD/2));
phi = linspace(phiF, alpha0, Nf);
length(edge1)
for i = 1: length(edge1)-1 % front from root
    xs = HubD*cosd(phi)/2;
    rs = HubD*sind(phi)/2;
    ys = linspace(edge1(3,i), edgeend(3,i), Nf);
    S1 = [xs; rs; ys];
    
    xs = HubD*cosd(phi)/2;
    rs = HubD*sind(phi)/2;
    ys = linspace(edge1(3,i+1), edgeend(3,i+1), Nf);
    S2 = [xs; rs; ys];
    
    S = [S1, S2];
    S(:,1:2:length(S)) = S1;
    S(:,2:2:length(S)) = S2;
    
    Vs = S';
    Fs = [1:length(S)-2; 2:length(S)-1; 3:length(S)]';
    
    F = [F; Fs+length(V)];
    V = [V; Vs];
    
    
    if i == 1
        STdw = [fliplr(S1(:,2:end)), STdw];
    elseif i == length(edge1)-1
        STup = [fliplr(S2(:,2:end)), STup];
    end
    
end

figure(88)
plot3(STup(1,:), STup(2,:), STup(3,:));
hold on;
plot3(STdw(1,:), STdw(2,:), STdw(3,:));
axis equal;

Nu = 15;
phi = linspace(alpha0, totalpha+alpha0, length(STup));
SSup = zeros(3,0); % store the shaft edge
for i = 1:length(STup) - 1 % hub upper face
    alpha = alpha0 + phi(i);
    xs = linspace(ShaftD*cosd(alpha)/2, STup(1,i), Nu);
    rs = linspace(ShaftD*sind(alpha)/2, STup(2,i), Nu);
    ys = 0*xs+STup(3,i);
    S1 = [xs; rs; ys];
    
    alpha = alpha0 + phi(i+1);
    xs = linspace(ShaftD*cosd(alpha)/2, STup(1,i+1), Nu);
    rs = linspace(ShaftD*sind(alpha)/2, STup(2,i+1), Nu);
    ys = 0*xs+STup(3,i+1);
    S2 = [xs; rs; ys];
    
    S = [S1, S2];
    S(:,1:2:length(S)) = S1;
    S(:,2:2:length(S)) = S2;
    
    Vs = S';
    Fs = [1:length(S)-2; 2:length(S)-1; 3:length(S)]';
    
    F = [F; Fs+length(V)];
    V = [V; Vs];
    
    SSup = [SSup, S1(:,1)];
    if i == length(STup)-1
        SSup = [SSup, S2(:,1)];
    end
    
end


phi = linspace(alpha0, totalpha+alpha0, length(STdw));
SSdw = zeros(3,0); % store the shaft edge
for i = 1:length(STdw) - 1 % hub down face
    alpha = alpha0 + phi(i);
    xs = linspace(ShaftD*cosd(alpha)/2, STdw(1,i), Nu);
    rs = linspace(ShaftD*sind(alpha)/2, STdw(2,i), Nu);
    ys = 0*xs+STdw(3,i);
    S1 = [xs; rs; ys];
    
    alpha = alpha0 + phi(i+1);
    xs = linspace(ShaftD*cosd(alpha)/2, STdw(1,i+1), Nu);
    rs = linspace(ShaftD*sind(alpha)/2, STdw(2,i+1), Nu);
    ys = 0*xs+STdw(3,i+1);
    S2 = [xs; rs; ys];
    
    S = [S1, S2];
    S(:,1:2:length(S)) = S1;
    S(:,2:2:length(S)) = S2;
    
    Vs = S';
    Fs = [1:length(S)-2; 2:length(S)-1; 3:length(S)]';
    
    F = [F; Fs+length(V)];
    V = [V; Vs];
    
    SSdw = [SSdw, S1(:,1)];
    if i == length(STdw)-1
        SSdw = [SSdw, S2(:,1)];
    end
end


Nh = 25;
nup = length(SSup);
ndw = length(SSdw);
diffn = abs(ndw-nup);
n = max(nup, ndw);
S2stor = NaN;
for i = 1:n-1 % the shaft inner face
    if mod(i,floor((n-2)/diffn)) == 0 % to fill in the missing points
        if nup > ndw
            SSdw = [SSdw(:,1:i), SSdw(:,i:end)]; % (the ith point repeats here)
        elseif ndw > nup
            SSup = [SSup(:,1:i), SSup(:,i:end)];
        end
    end
    
    ys = linspace(SSdw(3,i), SSup(3,i), Nh);
    xs = linspace(SSdw(1,i), SSup(1,i), Nh);
    rs = linspace(SSdw(2,i), SSup(2,i), Nh);
    S1 = [xs; rs; ys];
    
    ys = linspace(SSdw(3,i+1), SSup(3,i+1), Nh);
    xs = linspace(SSdw(1,i+1), SSup(1,i+1), Nh);
    rs = linspace(SSdw(2,i+1), SSup(2,i+1), Nh);
    S2 = [xs; rs; ys];
    
    if all(all(not(isnan(S2stor))))
        figure(88)
        plot3(S1(1,:), S1(2,:), S1(3,:), 'bo');
        hold on;
        plot3(S2stor(1,:), S2stor(2,:), S2stor(3,:), 'r*');
        S1 = S2stor;
        S2stor = NaN;
    end
    
    if mod(i,floor((n-2)/diffn)) == 0 
        if nup > ndw
            S1stor = S2;
            S1 = S1(:,floor(Nh/2):end);
            S2 = S2(:,floor(Nh/2)+1:end);
            S2stor = [S1stor(:,1:floor(Nh/2)), S2];
        elseif ndw > nup
            S1stor = S1;
            S1 = S1(:,1:floor(Nh/2)+1);
            S2 = S2(:,1:floor(Nh/2));
            S2stor = [S2, S1stor(:,floor(Nh/2)+1:end)];
        end
    end
    
    
    S = [S1, S2];
    if length(S1)>=length(S2)
        S(:,1:2:length(S)) = S1;
        S(:,2:2:length(S)) = S2;
    else
        S(:,1:2:length(S)) = S2;
        S(:,2:2:length(S)) = S1;
    end
    
    Vs = S';
    Fs = [1:length(S)-2; 2:length(S)-1; 3:length(S)]';
    
    F = [F; Fs+length(V)];
    V = [V; Vs];
end

figure(7)
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


% Remove duplicate vertices (the stl file will add them back again anyway though)
[V, indexm, indexn] =  unique(V, 'rows', 'stable');
F = indexn(F);

% Write to .stl CAD model
stlwrite('./Matlab_Out/test.stl', F, V, 'mode', 'ascii');

% mesh sizing analysis
figure(12)
plot(1:length(t), t, 'bo-')
hold on
plt = linspace(t(1), t(end), 200);
plot( linspace(1,length(t), 200), plt)
nint = zeros(length(t)-1, 1);
for i = 1: length(t)-1
    nint(i) = nnz(plt>t(i) & plt<t(i+1));
end
plot(2:length(t), nint);
grid on;

figure(7)
    
