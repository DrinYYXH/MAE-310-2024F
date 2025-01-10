clear;
clc;
 
% User-defined Parts
D_defined = 2;  % define the asumption for two-dimentional analysis plane stain(1) or plane stress(2)


% auto-going

% define Const.
E      = 1e9;
nu     = 0.3;
lambda = E*nu / ( (1 + nu) * (1 - 2*nu) );
miu    = E / (2*(1 + nu));

% stiffness matrix
if D_defined == 1
    D = [lambda + 2*miu , lambda         , 0;
         lambda         , lambda + 2*miu , 0;
         0              , 0              , miu];
    
elseif D_defined == 2
    D = [1              , nu             , 0;
         nu             , 1              , 0;
         0              , 0              , (1 - nu)/2] * (E / (1 - nu^2));
else
    error('Error: value of D_defined should be 1 or 2');
end


n_en   = 4;
n_el_x = 4;
n_el_y = 4;
n_sd   = 2;

g = @(x,y,i) 0;
h = @(x,y,i) 0;

f = @(x,y,i) (i == 1) * ((E*(nu - 4*x - 6*y - 2*nu*y + 4*x*y - 2*nu*x^2 + 2*x^2 + 4*y^2 + 4*nu*x*y + 1))/(2*(nu^2 - 1)) ) + ...
             (i == 2) * ((E*(nu - 6*x - 4*y - 2*nu*x + 4*x*y - 2*nu*y^2 + 4*x^2 + 2*y^2 + 4*nu*x*y + 1))/(2*(nu^2 - 1)) );

% exact solution
exact    = @(x,y,i) (i == 1) * (x*(1-x)*y*(1-y)) + (i == 2) * (x*(1-x)*y*(1-y));
exact_x  = @(x,y,i) (i == 1) * (1-2*x)*y*(1-y)   + (i == 2) * (1-2*x)*y*(1-y);
exact_y  = @(x,y,i) (i == 1) * x*(1-x)*(1-2*y)   + (i == 2) * x*(1-x)*(1-2*y);
exact_xx = @(x,y,i) (i == 1) * -2*y*(1-y)        + (i == 2) * -2*y*(1-y);
exact_yy = @(x,y,i) (i == 1) * -2*x*(1-x)        + (i == 2) * -2*x*(1-x);
exact_xy = @(x,y,i) (i == 1) * (1-2*x)*(1-2*y)   + (i == 2) * (1-2*x)*(1-2*y);
exact_yx = @(x,y,i) (i == 1) * (1-2*x)*(1-2*y)   + (i == 2) * (1-2*x)*(1-2*y);

% f = @(x,y,i) (i == 1) * ((2*E*y*(y - 1))/(nu^2 - 1) - (2*E*x*(nu/2 - 1/2)*(x - 1))/(nu^2 - 1) ) + ...
%              (i == 2) * ((E*(2*x - 1)*(2*y - 1))/(2*(nu - 1)) );
% 
% % exact solution
% exact    = @(x,y,i) (i == 1) * (x*(1-x)*y*(1-y)) + (i == 2) * 0;
% exact_x  = @(x,y,i) (i == 1) * (1-2*x)*y*(1-y)   + (i == 2) * 0;
% exact_y  = @(x,y,i) (i == 1) * x*(1-x)*(1-2*y)   + (i == 2) * 0;
% exact_xx = @(x,y,i) (i == 1) * -2*y*(1-y)        + (i == 2) * 0;
% exact_yy = @(x,y,i) (i == 1) * -2*x*(1-x)        + (i == 2) * 0;
% exact_xy = @(x,y,i) (i == 1) * (1-2*x)*(1-2*y)   + (i == 2) * 0;
% exact_yx = @(x,y,i) (i == 1) * (1-2*x)*(1-2*y)   + (i == 2) * 0;
 
u.exact    = exact;
u.exact_x  = exact_x;
u.exact_y  = exact_y;
u.exact_xx = exact_xx;
u.exact_yy = exact_yy;
u.exact_xy = exact_xy;
u.exact_yx = exact_yx;



mesh = Mesh4hole(n_en,n_sd);

% quadrature rule
n_int_xi  = 3;
n_int_eta = 3;
n_int     = n_int_xi * n_int_eta;
[xi, eta, weight] = Gauss2D(n_int_xi, n_int_eta);


displacement = FEM4hole(mesh,n_int,weight,xi,eta,f,D,g,h);


[e0,e1] = Error(mesh,displacement,n_int,weight,xi,eta,u);


%EOF