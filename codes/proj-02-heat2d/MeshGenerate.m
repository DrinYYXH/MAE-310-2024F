function mesh = MeshGenerate(n_en,n_el_x,n_el_y)
%MESHGENERATE 此处显示有关此函数的摘要
%   此处显示详细说明
% mesh generation
% n_en   = 3;               % number of nodes in an element
% n_el_x = 4;               % number of elements in x-dir
% n_el_y = 4;               % number of elements in y-dir
n_el   = n_el_x * n_el_y * 2; % total number of elements

n_np_x = n_el_x + 1;      % number of nodal points in x-dir
n_np_y = n_el_y + 1;      % number of nodal points in y-dir
n_np   = n_np_x * n_np_y; % total number of nodal points

x_coor = zeros(n_np, 1);
y_coor = x_coor;

hx = 1.0 / n_el_x;        % mesh size in x-dir
hy = 1.0 / n_el_y;        % mesh size in y-dir

% generate the nodal coordinates
for ny = 1 : n_np_y
  for nx = 1 : n_np_x
    index = (ny-1)*n_np_x + nx; % nodal index
    x_coor(index) = (nx-1) * hx;
    y_coor(index) = (ny-1) * hy;
  end
end

% IEN array
IEN = zeros(n_el, n_en);
for ex = 1 : n_el_x
  for ey = 1 : n_el_y
    ee = 2 * ((ey-1) * n_el_x + ex) - 1; % element index low triangle
    IEN(ee, 1) = (ey-1) * n_np_x + ex;
    IEN(ee, 2) = (ey-1) * n_np_x + ex + 1;
    IEN(ee, 3) =  ey    * n_np_x + ex;
    ee = 2 * ((ey-1) * n_el_x + ex); % element index up triangle
    IEN(ee, 1) = (ey-1) * n_np_x + ex + 1;
    IEN(ee, 2) =  ey    * n_np_x + ex + 1;
    IEN(ee, 3) =  ey    * n_np_x + ex;
  end
end

% ID array
ID = zeros(n_np,1);
counter = 0;
for ny = 2 : n_np_y - 1
  for nx = 2 : n_np_x - 1
    index = (ny-1)*n_np_x + nx;
    counter = counter + 1;
    ID(index) = counter;  
  end
end

n_eq = counter;

LM = ID(IEN);

mesh.IEN = IEN;
mesh.ID  = ID;
mesh.LM  = LM;

mesh.n_en   = n_en;
mesh.n_el   = n_el;
mesh.n_np   = n_np;
mesh.n_np_x = n_np_x;
mesh.n_np_y = n_np_y;
mesh.n_eq   = n_eq;
mesh.n_el_x = n_el_x;
mesh.n_el_y = n_el_y;

mesh.x_coor = x_coor;
mesh.y_coor = y_coor;

mesh.hx = hx;
mesh.hy = hy;

end

