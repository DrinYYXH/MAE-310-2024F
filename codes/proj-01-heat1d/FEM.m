function disp = FEM(mesh,h,g,n_int,f)
%FEM 此处显示有关此函数的摘要
%   此处显示详细说明

IEN    = mesh.IEN;
ID     = mesh.ID;
x_coor = mesh.coor;
hh     = mesh.hh;
n_el   = mesh.n_el;
n_en   = mesh.n_en;

pp = n_en - 1;
n_np = n_el * pp + 1;  % number of nodal points
n_eq = n_np - 1;       % number of equations

% allocate the stiffness matrix
K = zeros(n_eq, n_eq);
F = zeros(n_eq, 1);

% Setup the quadrature rule
[xi, weight] = Gauss(n_int, -1, 1);

for ee = 1 : n_el
  
  k_ele = zeros(n_en, n_en); % allocate a zero element stiffness matrix
  f_ele = zeros(n_en, 1);    % allocate a zero element load vector

  x_ele = x_coor(IEN(ee,:));

  for qua = 1 : n_int
    
    dx_dxi = 0.0;
    x_l = 0.0;
    for aa = 1 : n_en
      x_l = x_l + x_ele(aa) * PolyShape(pp, aa, xi(qua), 0);
      dx_dxi = dx_dxi + x_ele(aa) * PolyShape(pp, aa, xi(qua), 1);
    end
    dxi_dx = 1.0 / dx_dxi;

    for aa = 1 : n_en
      f_ele(aa) = f_ele(aa) + weight(qua) * PolyShape(pp, aa, xi(qua), 0) * f(x_l) * dx_dxi;
      for bb = 1 : n_en
        k_ele(aa, bb) = k_ele(aa, bb) + weight(qua) * PolyShape(pp, aa, xi(qua), 1) * PolyShape(pp, bb, xi(qua), 1) * dxi_dx;
      end
    end
  end
 
  % check the ID(IEN(ee, aa)) and ID(IEN(ee,bb)), if they are positive
  % put the element stiffness matrix into K
  for aa = 1 : n_en
    P = ID(IEN(ee,aa));
    if(P > 0)
      F(P) = F(P) + f_ele(aa);
      for bb = 1 : n_en
        Q = ID(IEN(ee,bb));
        if(Q > 0)
          K(P, Q) = K(P, Q) + k_ele(aa, bb);
        else
          F(P) = F(P) - k_ele(aa, bb) * g; % handles the Dirichlet boundary data
        end
      end
    end
  end

  if ee == 1
    F(ID(IEN(ee,1))) = F(ID(IEN(ee,1))) + h;
  end
end

% Solve Kd = F equation
d_temp = K \ F;

disp = [d_temp; g];


end

