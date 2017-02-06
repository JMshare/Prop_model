function [eg_from, eg_to, N] = blend_edges(eg_from, eg_to, part)


    n_from = length(eg_from);
    n_to = length(eg_to);
    N = max(n_from, n_to);
    n_diff = abs(n_from-n_to);
    if n_diff~=0
        sprintf('Blended edges in part %g.\n', part);
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
    

end