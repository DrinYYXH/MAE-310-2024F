clear;
clc;


syms r theta x y

Tx = 2;
R  = 1;



sigma_rr =  Tx/2 * (1 - (R^2)/(r^2)) + Tx/2 * (1 - 4 * (R^2)/(r^2) + 3 * (R^4)/(r^4)) * cos(2*theta);

sigma_tt =  Tx/2 * (1 + (R^2)/(r^2)) - Tx/2 * (1 + 3 * (R^4)/(r^4)) * cos(2*theta);

sigma_rt = -Tx/2 * (1 + 2 * (R^2)/(r^2) - 3 * (R^4)/(r^4)) * sin(2*theta);

T = [cos(theta)^2    , sin(theta)^2   , sin(theta)*cos(theta) ;
     sin(theta)^2    , cos(theta)^2   , -sin(theta)*cos(theta);
     -sin(2*theta)/2 , sin(2*theta)/2 , cos(2*theta)         ];


sigma = T * [sigma_rr ; sigma_tt ; sigma_rt];

r     = sqrt(x^2 + y^2);
theta = atan2(y, x);

f_x = -diff(sigma(1),x) - diff(sigma(3),y);
f_y = -diff(sigma(2),y) - diff(sigma(3),x);

simplify(f_x)
simplify(f_y)







%EOF
