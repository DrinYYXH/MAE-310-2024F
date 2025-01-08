clear;
clc;


syms x y E nu


 D = [1              , nu             , 0;
      nu             , 1              , 0;
      0              , 0              , (1 - nu)/2] * (E / (1 - nu^2));



exact_x  = x*(1-x)*y*(1-y);
exact_y  = x*(1-x)*y*(1-y);

exact_xx = (1-2*x)*y*(1-y);
exact_yy = (1-x)*x*(1-2*y);
exact_xy = exact_xx + exact_yy;


sigma = D * [exact_xx; exact_yy; exact_xy];

f_x = -diff(sigma(1),x) - diff(sigma(3),y);
f_y = -diff(sigma(2),y) - diff(sigma(3),x);

simplify(f_x)
simplify(f_y)

















 %EOF