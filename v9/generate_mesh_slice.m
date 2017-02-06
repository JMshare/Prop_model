function [Fs, Vs, F, V] = generate_mesh_slice(F, V, S1, S2)
    
    n1 = length(S1);
    n2 = length(S2);
    
    if n1 == n2
        S = [S1, S2];
        S(:,1:2:length(S)) = S1;
        S(:,2:2:length(S)) = S2;
    elseif n1 == n2+1
        S = [S1, S2];
        S(:,1:2:length(S)) = S1;
        S(:,2:2:length(S)) = S2;
    elseif n2 == n1+1
        S = [S2, S1];
        S(:,1:2:length(S)) = S2;
        S(:,2:2:length(S)) = S1;
    elseif n1 == n2+2
        S = [S1(:,1:end-1), S2];
        S(:,1:2:length(S)) = S1(:,1:end-1);
        S(:,2:2:length(S)) = S2;
        S = [S, S1(:,end)];
    elseif n2 == n1+2
        S = [S2(:,1:end-1), S1];
        S(:,1:2:length(S)) = S2(:,1:end-1);
        S(:,2:2:length(S)) = S1;
        S = [S, S2(:,end)];
    else
        sprintf('Error mesh slicing: given sices length differ by more than 2.');
        return;
    end
    
    Vs = S';
    Fs = [1:length(S)-2; 2:length(S)-1; 3:length(S)]';
    F = [F; Fs+length(V)];
    V = [V; Vs];

end