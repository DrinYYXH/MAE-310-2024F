clear all; clc;

kappa = 1.0; % conductivity

% exact solution
exact    = @(x,y) x*(1-x)*y*(1-y);
exact_x  = @(x,y) (1-2*x)*y*(1-y);
exact_y  = @(x,y) x*(1-x)*(1-2*y);
exact_xx = @(x,y) -2*y*(1-y);
exact_yy = @(x,y) -2*x*(1-x);
exact_xy = @(x,y) (1-2*x)*(1-2*y);
exact_yx = @(x,y) (1-2*x)*(1-2*y);

u.exact    = exact;
u.exact_x  = exact_x;
u.exact_y  = exact_y;
u.exact_xx = exact_xx;
u.exact_yy = exact_yy;
u.exact_xy = exact_xy;
u.exact_yx = exact_yx;


f = @(x,y) 2.0*kappa*x*(1-x) + 2.0*kappa*y*(1-y); % source term

% quadrature rule
% n_int_xi  = 3;
% n_int_eta = 3;
% n_int     = n_int_xi * n_int_eta;
% [xi, eta, weight] = Gauss2D(n_int_xi, n_int_eta);

% quadrature rule for triangle
n_int = 3;
weight = [1/3 1/3 1/3];
xi     = [1/6 2/3 1/6];
eta    = [1/6 1/6 2/3];

%% mesh generation
n_en   = 3;               % number of nodes in an element
n_el_x = 60;               % number of elements in x-dir
n_el_y = 60;               % number of elements in y-dir
% n_el   = n_el_x * n_el_y * 2; % total number of elements
% 
% n_np_x = n_el_x + 1;      % number of nodal points in x-dir
% n_np_y = n_el_y + 1;      % number of nodal points in y-dir
% n_np   = n_np_x * n_np_y; % total number of nodal points
% 
% x_coor = zeros(n_np, 1);
% y_coor = x_coor;
% 
% hx = 1.0 / n_el_x;        % mesh size in x-dir
% hy = 1.0 / n_el_y;        % mesh size in y-dir


mesh = MeshGenerate(n_en,n_el_x,n_el_y);

IEN = mesh.IEN;
ID  = mesh.ID;
LM  = mesh.LM;

n_el   = mesh.n_el;
n_np   = mesh.n_np;
n_np_x = mesh.n_np_x;
n_np_y = mesh.n_np_y;
n_eq   = mesh.n_eq;

x_coor = mesh.x_coor;
y_coor = mesh.y_coor;

hx = mesh.hx;
hy = mesh.hy;


%% FEM

displacement = FEM(mesh,n_int,weight,xi,eta,f,kappa);




%% error

[e0,e1] = Error(mesh,displacement,n_int,weight,xi,eta,u);




%% error visualization
clear;
clc;

n_en   = 3.0;                                  % number of nodes in an element

n_ele   = [2,4,6,8,10,12,14,16,18,20];
n_ele_x = n_ele;               % number of elements in x-dir
n_ele_y = n_ele;               % number of elements in y-dir


% quadrature rule for triangle
n_int = 3;
weight = [1/3 1/3 1/3];
xi     = [1/6 2/3 1/6];
eta    = [1/6 1/6 2/3];

e0 = zeros(1,length(n_ele));
e1 = zeros(1,length(n_ele));


for i = 1 : length(n_ele)

kappa = 1.0; % conductivity

    % exact solution
    exact    = @(x,y) x*(1-x)*y*(1-y);
    exact_x  = @(x,y) (1-2*x)*y*(1-y);
    exact_y  = @(x,y) x*(1-x)*(1-2*y);
    exact_xx = @(x,y) -2*y*(1-y);
    exact_yy = @(x,y) -2*x*(1-x);
    exact_xy = @(x,y) (1-2*x)*(1-2*y);
    exact_yx = @(x,y) (1-2*x)*(1-2*y);
    
    u.exact    = exact;
    u.exact_x  = exact_x;
    u.exact_y  = exact_y;
    u.exact_xx = exact_xx;
    u.exact_yy = exact_yy;
    u.exact_xy = exact_xy;
    u.exact_yx = exact_yx;
    
    f = @(x,y) 2.0*kappa*x*(1-x) + 2.0*kappa*y*(1-y); % source term

    % Setup the mesh

    n_el_x = n_ele_x(i); 
    n_el_y = n_ele_y(i); 

    %initial the mesh for equally spaced nodes
    mesh = MeshGenerate(n_en,n_el_x,n_el_y);

    IEN = mesh.IEN;
    ID  = mesh.ID;
    LM  = mesh.LM;
    
    n_el   = mesh.n_el;
    n_np   = mesh.n_np;
    n_np_x = mesh.n_np_x;
    n_np_y = mesh.n_np_y;
    n_eq   = mesh.n_eq;
    
    x_coor = mesh.x_coor;
    y_coor = mesh.y_coor;
    
    hx = mesh.hx;
    hy = mesh.hy;

    %solve FEM problem
    displacement = FEM(mesh,n_int,weight,xi,eta,f,kappa);


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












% EOF