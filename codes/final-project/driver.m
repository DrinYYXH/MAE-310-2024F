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



mesh = MeshGenerate(n_en,n_el_x,n_el_y,n_sd);

% quadrature rule
n_int_xi  = 3;
n_int_eta = 3;
n_int     = n_int_xi * n_int_eta;
[xi, eta, weight] = Gauss2D(n_int_xi, n_int_eta);


displacement = FEM(mesh,n_int,weight,xi,eta,f,D,g,h);


[e0,e1] = Error(mesh,displacement,n_int,weight,xi,eta,u);


%% error visualization
clear;
clc;

D_defined = 2;

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
n_sd   = 2;

n_ele   = [8,10,12,14,16,18,20];
n_ele_x = n_ele;               % number of elements in x-dir
n_ele_y = n_ele;               % number of elements in y-dir




% quadrature rule
n_int_xi  = 3;
n_int_eta = 3;
n_int     = n_int_xi * n_int_eta;
[xi, eta, weight] = Gauss2D(n_int_xi, n_int_eta);



f = @(x,y,i) (i == 1) * ((E*(nu - 4*x - 6*y - 2*nu*y + 4*x*y - 2*nu*x^2 + 2*x^2 + 4*y^2 + 4*nu*x*y + 1))/(2*(nu^2 - 1)) ) + ...
             (i == 2) * ((E*(nu - 4*x - 6*y - 2*nu*y + 4*x*y - 2*nu*x^2 + 2*x^2 + 4*y^2 + 4*nu*x*y + 1))/(2*(nu^2 - 1)) );


% exact solution
exact    = @(x,y,i) (i == 1) * (x*(1-x)*y*(1-y)) + (i == 2) * (x*(1-x)*y*(1-y));
exact_x  = @(x,y,i) (i == 1) * (1-2*x)*y*(1-y)   + (i == 2) * (1-2*x)*y*(1-y);
exact_y  = @(x,y,i) (i == 1) * x*(1-x)*(1-2*y)   + (i == 2) * x*(1-x)*(1-2*y);
exact_xx = @(x,y,i) (i == 1) * -2*y*(1-y)        + (i == 2) * -2*y*(1-y);
exact_yy = @(x,y,i) (i == 1) * -2*x*(1-x)        + (i == 2) * -2*x*(1-x);
exact_xy = @(x,y,i) (i == 1) * (1-2*x)*(1-2*y)   + (i == 2) * (1-2*x)*(1-2*y);
exact_yx = @(x,y,i) (i == 1) * (1-2*x)*(1-2*y)   + (i == 2) * (1-2*x)*(1-2*y);

u.exact    = exact;
u.exact_x  = exact_x;
u.exact_y  = exact_y;
u.exact_xx = exact_xx;
u.exact_yy = exact_yy;
u.exact_xy = exact_xy;
u.exact_yx = exact_yx;

e0 = zeros(1,length(n_ele));
e1 = zeros(1,length(n_ele));


for i = 1 : length(n_ele)


    % Setup the mesh

    n_el_x = n_ele_x(i); 
    n_el_y = n_ele_y(i); 

    %initial the mesh for equally spaced nodes
    mesh = MeshGenerate(n_en,n_el_x,n_el_y,n_sd);

    %unpack mesh
    IEN = mesh.IEN;
    ID  = mesh.ID;
    LM  = mesh.LM;
    
    n_en   = mesh.n_en;
    n_el   = mesh.n_el;
    n_np   = mesh.n_np;
    n_np_x = mesh.n_np_x;
    n_np_y = mesh.n_np_y;
    n_eq   = mesh.n_eq;
    n_el_x = mesh.n_el_x;
    n_el_y = mesh.n_el_y;
    n_sd   = mesh.n_sd;
    
    n_ee   = n_en * n_sd;
    
    x_coor = mesh.x_coor;
    y_coor = mesh.y_coor;
    
    hx = mesh.hx;
    hy = mesh.hy;

    %solve FEM problem
    displacement = FEM(mesh,n_int,weight,xi,eta,f,D);


    [e0(i),e1(i)] = Error(mesh,displacement,n_int,weight,xi,eta,u);


end


hh = (1./n_ele)*sqrt(2.0);

figure(1)
plot(log(hh),log(e0),'-xr');
hold on;
plot(log(hh),log(e1),'-.k');

xlabel('log(h)');
ylabel('log(error)');
title('log-log h-error');
legend('e_{0}', 'e_{1}','Location','best');

p1 = polyfit(log(hh),log(e0),1);
p2 = polyfit(log(hh),log(e1),1);


slope1 = sprintf('Slope e_{0} = %.2f', p1(1));
slope2 = sprintf('Slope e_{1} = %.2f', p2(1));
% 根据需要调整文字位置
text(log(hh(end))+0.5, log(e0(end)), slope1, 'FontSize', 12, 'Color', 'r');
text(log(hh(end))+0.5, log(e1(end)), slope2, 'FontSize', 12, 'Color', 'k');

hold off;










%EOF