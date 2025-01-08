clear;
clc;


syms x y E nu


 D = [1              , nu             , 0
      nu             , 1              , 0
      0              , 0              , (1 - nu)/2] * (E / (1 - nu^2));



exact_x  = x*(x-1)*y*(y-1);
exact_y  = x*(x-1)*y*(y-1);

exact_xx = (2*x - 1)*y*(y-1);
exact_yy = (2*y - 1)*x*(x-1);
exact_xy = exact_xx + exact_yy;


sigma = D * [exact_xx; exact_yy; exact_xy];

f_x = -diff(sigma(1),x) - diff(sigma(3),y);
f_y = -diff(sigma(2),y) - diff(sigma(3),x);



















 %EOF