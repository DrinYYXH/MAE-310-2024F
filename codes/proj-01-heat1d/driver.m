clear all; clc; %clf; % clean the memory, screen, and figure

% Problem definition
f = @(x) -20*x.^3; % f(x) is the source
g = 1.0;           % u    = g  at x = 1
h = 0.0;           % -u,x = h  at x = 0
u = @(x) x.^5;      % test solve
u_x = @(x) 5*x.^4;

% Setup the mesh
pp   = 1;              % polynomial degree
n_en = pp + 1;         % number of element or local nodes
n_el = 2;              % number of elements
n_np = n_el * pp + 1;  % number of nodal points
n_eq = n_np - 1;       % number of equations
n_int = 10;

hh = 1.0 / (n_np - 1); % space between two adjacent nodes
x_coor = 0 : hh : 1;   % nodal coordinates for equally spaced nodes

[ID,IEN] = GenerateMesh(n_el,n_en);

% Setup the quadrature rule
[xi, weight] = Gauss(n_int, -1, 1);

% allocate the stiffness matrix
K = spalloc(n_eq, n_eq, (2*pp+1)*n_eq);
F = zeros(n_eq, 1);

% Assembly of the stiffness matrix and load vector
for ee = 1 : n_el
  k_ele = zeros(n_en, n_en); % allocate a zero element stiffness matrix
  f_ele = zeros(n_en, 1);    % allocate a zero element load vector

  x_ele = x_coor(IEN(ee,:)); % x_ele(aa) = x_coor(A) with A = IEN(aa, ee)

  % quadrature loop
  for qua = 1 : n_int    
    dx_dxi = 0.0;
    x_l = 0.0;
    for aa = 1 : n_en
      x_l    = x_l    + x_ele(aa) * PolyShape(pp, aa, xi(qua), 0);
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
 
  % Assembly of the matrix and vector based on the ID or LM data
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
end

% ee = 1 F = NA(0)xh
F(ID(IEN(1,1))) = F(ID(IEN(1,1))) + h;

% Solve Kd = F equation
d_temp = K \ F;

disp = [d_temp; g];

% Postprocessing: visualization
%plot(x_coor, disp, '--r','LineWidth',3);

%x_sam = 0 : 0.01 : 1;
%y_sam = x_sam.^5;
%hold on;
%plot(x_sam, y_sam, '-k', 'LineWidth', 3);


%% visualization prof


n_sam = 20;
xi_sam = -1 : (2/n_sam) : 1;

x_sam = zeros(n_el * n_sam + 1, 1);
y_sam = x_sam; % store the exact solution value at sampling points
u_sam = x_sam; % store the numerical solution value at sampling pts

for ee = 1 : n_el
  x_ele = x_coor( IEN(ee, :) );
  u_ele = disp( IEN(ee, :) );

  if ee == n_el
    n_sam_end = n_sam+1;
  else
    n_sam_end = n_sam;
  end

  for ll = 1 : n_sam_end
    x_l = 0.0;
    u_l = 0.0;
    for aa = 1 : n_en
      x_l = x_l + x_ele(aa) * PolyShape(pp, aa, xi_sam(ll), 0);
      u_l = u_l + u_ele(aa) * PolyShape(pp, aa, xi_sam(ll), 0);
    end

    x_sam( (ee-1)*n_sam + ll ) = x_l;
    u_sam( (ee-1)*n_sam + ll ) = u_l;
    y_sam( (ee-1)*n_sam + ll ) = x_l^5;
  end
end


plot(x_sam, u_sam, '-r','LineWidth',3);
hold on;
plot(x_sam, y_sam, '-k','LineWidth',3);





%% visualization lzy

% maybe a stupid way that make a large array just call it IP array, about the same as IEN array
% IP array is mainly the find every point of uu from partial coor to total coor

tt = 5;                                 %divide the points in each partial element
ipn = 1 : 1/tt : n_en;                  %make each points a mark

dd = (tt - 1)*(n_en - 1) + n_en;        %in order to make IP easy
IP = zeros(n_en,length(ipn));
for ee = 1 : n_el
    for xip = 1 : length(ipn)
        IP(ee,xip) = xip + (dd - 1)*(ee - 1);
    end
end


uu = zeros(IP(n_el,length(ipn)),1);     %u_h actually
x_xii = zeros(IP(n_el,length(ipn)),1);  %the x of each u_h in total coor

for ee = 1 : n_el
    for xip = 1 : length(ipn)
        for aa = 1 : n_en
        xii = 2*(xip - 1)/(length(ipn) - 1) - 1;

        x_xii(IP(ee,xip)) = x_xii(IP(ee,xip)) + x_coor(IEN(ee,aa)) * PolyShape(pp,aa,xii,0);

        uu(IP(ee,xip)) = uu(IP(ee,xip)) + disp(IEN(ee,aa)) * PolyShape(pp,aa,xii,0);

        end

        if (ee > 1) && (xip == 1)
            x_xii(IP(ee,xip)) = x_xii(IP(ee,xip))/2;
            uu(IP(ee,xip)) = uu(IP(ee,xip))/2;
        end

    end
end

% draw a picture of u_h-x compare with u-x
figure(1)
plot(x_xii,uu,'xr');
hold on;
plot(x_xii,x_xii.^5,'k');


%% error H1 L2

[e_L2,e_H1] = Error(n_el,n_en,IEN,disp,n_int,u,u_x);
















% EOF