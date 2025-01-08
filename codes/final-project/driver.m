clear;
clc;
 
%% User-defined Parts
D_defined = 2;  % define the asumption for two-dimentional analysis plane stain(1) or plane stress(2)








%% auto-going

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
n_el_x = 10;
n_el_y = 10;
n_sd   = 2;


exact = @(x,y,i) (i == 1) * (x*(x-1)*y*(y-1)) + (i == 2) * (x*(x-1)*y*(y-1));

f = @(x,y,i) (i == 1) * ((E*(nu - 4*x - 6*y - 2*nu*y + 4*x*y - 2*nu*x^2 + 2*x^2 + 4*y^2 + 4*nu*x*y + 1))/(2*(nu^2 - 1)) ) + ...
             (i == 2) * ((E*(nu - 4*x - 6*y - 2*nu*y + 4*x*y - 2*nu*x^2 + 2*x^2 + 4*y^2 + 4*nu*x*y + 1))/(2*(nu^2 - 1)) );




mesh = MeshGenerate(n_en,n_el_x,n_el_y,n_sd);

% quadrature rule
n_int_xi  = 3;
n_int_eta = 3;
n_int     = n_int_xi * n_int_eta;
[xi, eta, weight] = Gauss2D(n_int_xi, n_int_eta);


displacement = FEM(mesh,n_int,weight,xi,eta,f,D);



















%EOF