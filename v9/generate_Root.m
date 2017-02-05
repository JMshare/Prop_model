function [r, c, beta, rake, skew, thickness, A] = generate_Root(r, c, beta, rake, skew, thickness, A, dRoot, rRoot, cRoot, tRoot, betaRoot, skewRoot, thickRoot, rakeRoot, Nrt)

    incr = dRoot/Nrt;
    Y = cell2mat(A(1));
    Y = Y(:,2);
    tstart = abs(Y);
    Ydd = cell2mat(A(2));
    Ydd = Ydd(:,2);
    tstartt = abs(Ydd);
    Yddd = cell2mat(A(3));
    Yddd = Yddd(:,2);
    tstarttt = abs(Yddd);
    r_Root = [];
    c_Root = [];
    beta_Root = [];
    rake_Root = [];
    skew_Root = [];
    thick_Root = [];
    A_Root = A(1);
    drstart = -(r(1)-r(2));
    for i = 1:Nrt
        r_Root = [r_Root; Lagrn3(i*incr, [0,r(1)], [dRoot, rRoot], [-drstart, r(2)], [-2*drstart, r(3)])]; % 
        v = Lagrn1(i*incr, [0, 2*r(1)*sin(c(1)/(r(1)*2))], [dRoot, 2*rRoot*sin(cRoot/(rRoot*2))]);
        c_Root = [c_Root; 2*r_Root(end)*asin(v/(2*r_Root(end)))]; %  
        beta_Root = [beta_Root; Lagrn1(i*incr, [0, beta(1)], [dRoot, betaRoot])];
        rake_Root = [rake_Root; Lagrn1(i*incr, [0, rake(1)], [dRoot, rakeRoot])];
        skew_Root = [skew_Root; Lagrn1(i*incr, [0, skew(1)], [dRoot+rRoot, skewRoot])];
        thick_Root = [thick_Root; Lagrn1(i*incr, [0, thickness(1)], [dRoot, thickRoot])];
        
        Y = cell2mat(A_Root(end));
        X = Y(:,1);
        Y = Y(:,2);
        for j = (ceil(length(Y)/2)+i) : (length(Y)-i)
            t = Lagrn1(i*incr, [0, tstart(j)], [dRoot, tRoot]); %
            if (Y(j) < 0) 
                Y(j) = -t;
            end
        end
        A_Root(i) = {[X, Y]};
    end
    r = [flipud(r_Root); r];
    c = [flipud(c_Root); c];
    beta = [flipud(beta_Root); beta];
    rake = [flipud(rake_Root); rake];
    skew = [flipud(skew_Root); skew];
    thickness = [flipud(thick_Root); thickness];
    A = [flipud(A_Root'); A];


end