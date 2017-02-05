function [F, V, eg_st, eg_nd] = mesh_surface(Xs, Rs, Ys, Ns, F, V, close, plti, fign1, fign2)
    
    if close == 1
        Xs = [Xs, Xs(:,1)];
        Rs = [Rs, Rs(:,1)];
        Ys = [Ys, Ys(:,1)];
    end
    
    N = size(Xs, 2);
    Fa = zeros(0,3);
    Va = zeros(0,3);
    for i = 1:N-1
        [S1] = spline_cubic3([Xs(:,i)'; Rs(:,i)'; Ys(:,i)'], Ns);
        adapt = 0;
        %adapt = 1*(eg_to(3,i+1)>eg_to(3,i)) - 1*(eg_to(3,i+1)<eg_to(3,i));
        [S2] = spline_cubic3([Xs(:,i+1)'; Rs(:,i+1)'; Ys(:,i+1)'], Ns+adapt);

        [Fs, Vs, Fa, Va] = generate_mesh_slice(Fa, Va, S1, S2);

        Ns = length(S2);
        
        if i == 1
            eg_st = S1;
        end
        if i == N-1
            eg_nd = S2;
        end
        
        if i == plti
            figure(fign1)
            axis equal;
            hold on;
            trimesh(Fs, Vs(:,1),Vs(:,2), Vs(:,3));
            title('The whole observed slice');
            figure(fign2)
            trimesh(Fs(end-3:end,:), Vs(:,1),Vs(:,2), Vs(:,3), 'FaceColor', 'blue');
            axis equal;
            hold on;
        end
    
    end
    
    F = [F; Fa+length(V)];
    V = [V; Va];
    
end