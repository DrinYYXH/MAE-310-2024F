clear all; clc; % clean the memory and screen

% Define the external source or force and boundary data
f = @(x) -20*x.^3; % f(x) = x
g = 1.0;           % u    = g  at x = 1
h = 0.0;           % -u,x = h  at x = 0

% Setup the mesh
pp   = 1;              % polynomial degree
n_en = pp + 1;         % number of element or local nodes
n_el = 4;              % number of elements
n_np = n_el * pp + 1;  % number of nodal points
n_eq = n_np - 1;       % number of equations
n_int = 10;

hh = 1.0 / (n_np - 1); % space between two adjacent nodes
x_coor = 0 : hh : 1;   % nodal coordinates for equally spaced nodes

IEN = zeros(n_el, n_en);

for ee = 1 : n_el
  for aa = 1 : n_en
    IEN(ee, aa) = (ee - 1) * pp + aa;
  end
end

% Setup the ID array for the problem
ID = 1 : n_np;
ID(end) = 0;

% Setup the quadrature rule
[xi, weight] = Gauss(n_int, -1, 1);

% allocate the stiffness matrix
K = zeros(n_eq, n_eq);
F = zeros(n_eq, 1);

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


%% 

% X = 0 : hh : 1;
% uu = zeros(size(X));
% Xi = 2*(mod(X,hh)/hh) - 1;
% 
% for xx = 1 : length(X)
%     for aa = 1 : n_en
%         p_temp = ceil(X(xx)/hh);
%         if p_temp < 1
%             ee = 1;
%         else
%             ee = ceil(X(xx)/hh);
%         end
% 
%         % syms xii
%         % x_xii = x_coor(IEN(ee,aa))*PolyShape(pp,aa,xii,0);
%         % xii_x = finverse(x_xii,xii);
% 
%         uu(xx) = uu(xx) + disp(IEN(ee,aa)) * PolyShape(pp,aa,Xi(xx),0);
%     end
% end
% 
% 
% plot(X,uu);

%%

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








% EOF