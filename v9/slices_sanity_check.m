function [S1, S2] = slices_sanity_check(S1, S2)


    if all(S2(:,end) == S1(:,end)) % if blended, you want one less element there.
        S2 = S2(:,1:end-1);
    end
    if all(S2(:,1) == S1(:,1)) % if blended, you want one less element there.
        S2 = S2(:,2:end);
    end

end