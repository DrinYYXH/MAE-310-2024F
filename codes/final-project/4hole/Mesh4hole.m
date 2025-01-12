function mesh = Mesh4hole(n_en,n_sd)
%MESHGENERATE 此处显示有关此函数的摘要
%   此处显示详细说明

addpath('H:\CSM\Code_HW\MAE-310-2024F\codes\final-project');

%% data input
% gmesh progress
% run('gmesh_output.m');
% run('gmesh_input.m');
run('output_hole.m');

n_np   = msh.nbNod - 1;

x_coor = zeros(n_np, 1);
y_coor = x_coor;

quads    = msh.QUADS;
lines    = msh.LINES;
position = msh.POS;

n_el = length(quads);

%% IEN array
IEN = zeros(n_el, n_en);
for ee = 1 : n_el
    for aa = 1 : n_en
    IEN(ee,aa) = quads(ee,aa);
    end
end



%% ID array
ID = zeros(n_np,n_sd);

BC = [1, 2, 3, 4;
      1, 2, 3, 4];

ID_abandon = zeros(n_np,n_sd);
Dirichlet_BC = [1, 2, 3, 4;
                1, 2, 3, 4];
for i = 1 : length(lines)
    for j = 1 : n_sd
    is_present = ismember(lines(i,3),Dirichlet_BC(j,:));
        if is_present
            ID_abandon(lines(i,1),j) = 1;
        end
    end
end

mesh.ID_abandon = ID_abandon;

counter = 0;
for i = 1 : n_np
    for j = 1 : n_sd
        if ID_abandon(i,j) == 0
            counter = counter + 1;
            ID(i,j) = counter;
        end
    end
end

n_eq = counter;



%% generate the nodal coordinates
for ee = 1 : n_el
  for aa = 1 : n_en
    x_coor(IEN(ee,aa)) = position(IEN(ee,aa),1);
    y_coor(IEN(ee,aa)) = position(IEN(ee,aa),2);
  end
end



%% IDH for Neumann B.C.
Neumann_BC = BC - Dirichlet_BC;
IDH = zeros(n_np,n_sd);
for i = 1 : length(lines)
    for j = 1 : n_sd
    is_present = ismember(lines(i,3),Neumann_BC(j,:));
        if is_present
            IDH(lines(i,1),j) = 1;
        end
    end
end


%% pack mesh
mesh.IDH  = IDH;

LM = ID(IEN);

mesh.IEN = IEN;
mesh.ID  = ID;
mesh.LM  = LM;



mesh.n_en   = n_en;
mesh.n_el   = n_el;
mesh.n_np   = n_np;
mesh.n_eq   = n_eq;
mesh.n_sd   = n_sd;

mesh.x_coor = x_coor;
mesh.y_coor = y_coor;

%EOF


