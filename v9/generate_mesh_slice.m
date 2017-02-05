function [Fs, Vs, F, V] = generate_mesh_slice(F, V, S1, S2)
    

    if length(S1) >= length(S2)
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
    F = [F; Fs+length(V)];
    V = [V; Vs];

end