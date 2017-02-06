function [F, V, eg_st, eg_nd] = mesh_area(eg_from, eg_to, Ns, F, V)
    
    N = length(eg_from);
    Fa = zeros(0,3);
    Va = zeros(0,3);
    for i = 1:N-1 
        S1 = spline_cylinder3([eg_from(:,i)'; eg_to(:,i)'],...
                         Ns);
        adapt = 0;
        %adapt = 1*(eg_to(3,i+1)>eg_to(3,i)) - 1*(eg_to(3,i+1)<eg_to(3,i));
        Ns = Ns+adapt;
        S2 = spline_cylinder3([eg_from(:,i+1)'; eg_to(:,i+1)'],... 
                         Ns);
        
        [S1, S2] = slices_sanity_check(S1, S2);

        [~, ~, Fa, Va] = generate_mesh_slice(Fa, Va, S1, S2);

        if i == 1
            eg_st = S1;
        end
        if i == N-1
            eg_nd = S2;
        end
    end
    
    F = [F; Fa+length(V)];
    V = [V; Va];
    
end