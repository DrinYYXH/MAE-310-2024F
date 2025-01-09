function mesh = MeshGenerate(n_en,n_el_x,n_el_y,n_sd)
%MESHGENERATE 此处显示有关此函数的摘要
%   此处显示详细说明

% gmesh progress
% run('gmesh_output.m');
% run('gmesh_input.m');
run('output.m');


n_el   = n_el_x * n_el_y;

n_np_x = n_el_x + 1;      % number of nodal points in x-dir
n_np_y = n_el_y + 1;      % number of nodal points in y-dir
n_np   = n_np_x * n_np_y; % total number of nodal points

x_coor = zeros(n_np, 1);
y_coor = x_coor;

hx = 1.0 / n_el_x;        % mesh size in x-dir
hy = 1.0 / n_el_y;        % mesh size in y-dir

quads    = msh.QUADS;
lines    = msh.LINES;
position = msh.POS;
length(lines);


% IEN array
IEN = zeros(n_el, n_en);
for ee = 1 : n_el
    for aa = 1 : n_en
    IEN(ee,aa) = quads(ee,aa);
    end
end

% ID array
ID = zeros(n_np,n_sd);

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


% generate the nodal coordinates
for ee = 1 : n_el
  for aa = 1 : n_en
    x_coor(IEN(ee,aa)) = position(IEN(ee,aa),1);
    y_coor(IEN(ee,aa)) = position(IEN(ee,aa),2);
  end
end


mesh.ID_abandon = ID_abandon;




% % total number of elements
% if n_en == 3
%     n_el   = n_el_x * n_el_y * 2; 
%  elseif n_en == 4
%     n_el   = n_el_x * n_el_y;
% else
%     error('Error: value of a should be 3 or 4.');
% end     
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
% 
% % generate the nodal coordinates
% for ny = 1 : n_np_y
%   for nx = 1 : n_np_x
%     index = (ny-1)*n_np_x + nx; % nodal index
%     x_coor(index) = (nx-1) * hx;
%     y_coor(index) = (ny-1) * hy;
%   end
% end
% 
% % IEN array
% IEN = zeros(n_el, n_en);
% 
% if n_en == 3
%     for ex = 1 : n_el_x
%       for ey = 1 : n_el_y
%         ee = 2 * ((ey-1) * n_el_x + ex) - 1; % element index low triangle
%         IEN(ee, 1) = (ey-1) * n_np_x + ex;
%         IEN(ee, 2) = (ey-1) * n_np_x + ex + 1;
%         IEN(ee, 3) =  ey    * n_np_x + ex;
%         ee = 2 * ((ey-1) * n_el_x + ex); % element index up triangle
%         IEN(ee, 1) = (ey-1) * n_np_x + ex + 1;
%         IEN(ee, 2) =  ey    * n_np_x + ex + 1;
%         IEN(ee, 3) =  ey    * n_np_x + ex;
%       end
%     end
% elseif n_en == 4
%      for ex = 1 : n_el_x
%       for ey = 1 : n_el_y
%         ee = (ey-1) * n_el_x + ex; % element index
%         IEN(ee, 1) = (ey-1) * n_np_x + ex;
%         IEN(ee, 2) = (ey-1) * n_np_x + ex + 1;
%         IEN(ee, 3) =  ey    * n_np_x + ex + 1;
%         IEN(ee, 4) =  ey    * n_np_x + ex;
%       end
%      end
% else
%     error('Error: value of a should be 3 or 4.');
% end
% 
% % ID array
% ID = zeros(n_np,n_sd);
% counter = 0;
% for ny = 2 : n_np_y - 1
%   for nx = 2 : n_np_x - 1
%     index = (ny-1)*n_np_x + nx;
%     for i = 1 : n_sd
%         counter = counter + 1;
%         ID(index,i) = counter;
%     end
%   end
% end
% 
% n_eq = counter;


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
mesh.n_sd   = n_sd;

mesh.x_coor = x_coor;
mesh.y_coor = y_coor;

mesh.hx = hx;
mesh.hy = hy;

% end

