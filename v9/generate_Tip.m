function [r, c, beta, rake, skew, thickness, A] = generate_Tip(r, c, beta, rake, skew, thickness, A, dTip, rTip, cTip, tTip, skewTip, Nt)

    incr = dTip/Nt;
    Y = cell2mat(A(end));
    Y = Y(:,2);
    tend = abs(Y);
    Ydd = cell2mat(A(end-1));
    Ydd = Ydd(:,2);
    tendd = abs(Ydd);
    Yddd = cell2mat(A(end-2));
    Yddd = Yddd(:,2);
    tenddd = abs(Yddd);
    r_Tip = [];
    c_Tip = [];
    beta_Tip = [];
    rake_Tip = [];
    skew_Tip = [];
    thick_Tip = [];
    A_Tip = A(end);
    drend = r(end) - (r(end-1) - dTip);
    for i = 1:Nt
        r_Tip = [r_Tip; Lagrn3(i*incr, [0,r(end)], [dTip, rTip], [2*dTip, r(end)], [-drend, r(end-1)])]; % quadrtic mesh densing to match the quadratic c reduction
        c_Tip = [c_Tip; Lagrn3(i*incr, [0, c(end)], [dTip, cTip], [-drend, c(end-1)], [-2*drend, c(end-2)])]; % quadratic chord reduction
        beta_Tip = [beta_Tip; Lagrn3(i*incr, [0, beta(end)], [-drend, beta(end-1)], [-2*drend, beta(end-2)], [-3*drend, beta(end-3)])];
        rake_Tip = [rake_Tip; Lagrn3(i*incr, [0, rake(end)], [-drend, rake(end-1)], [-2*drend, rake(end-2)], [-3*drend, rake(end-3)])];
        skew_Tip = [skew_Tip; Lagrn3(i*incr, [0, skew(end)], [dTip, skewTip], [-drend, skew(end-1)], [-2*drend, skew(end-2)])];
        thick_Tip = [thick_Tip; Lagrn3(i*incr, [0, thickness(end)], [-drend, thickness(end-1)], [-2*drend, thickness(end-2)], [-3*drend, thickness(end-3)])];
        
        Y = cell2mat(A_Tip(end));
        X = Y(:,1);
        Y = Y(:,2);
        for j = 1:length(Y)
            t = Lagrn3(i*incr, [0, tend(j)], [dTip, tTip], [-drend, tendd(j)], [-2*drend, tenddd(j)]); % quadratic airfoil thickness reduction
            if abs(Y(j)) > t
                Y(j) = sign(Y(j))*t;
            end
        end
        A_Tip(i) = {[X, Y]};
    end
    r = [r; r_Tip];
    c = [c; c_Tip];
    beta = [beta; beta_Tip];
    rake = [rake; rake_Tip];
    skew = [skew; skew_Tip];
    thickness = [thickness; thick_Tip];
    A = [A; A_Tip'];


end