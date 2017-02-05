function [y] = Lagrn1(x, p1, p2)

x1 = p1(1);
y1 = p1(2);
x2 = p2(1);
y2 = p2(2);


y = (y1*(x-x2))/((x1-x2)) +...
    (y2*(x-x1))/((x2-x1));

end