function [e_L2,e_H1] = Error(mesh,disp,n_int,u,u_x)
%ERROR 此处显示有关此函数的摘要
%   此处显示详细说明

IEN    = mesh.IEN;
ID     = mesh.ID;
x_coor = mesh.coor;
hh     = mesh.hh; 
n_el   = mesh.n_el;
n_en   = mesh.n_en;

pp = n_en - 1;         % polynomial degree

% Setup the quadrature rule
[xi, weight] = Gauss(n_int, -1, 1);


%e_L2
e_L2_up   = 0;
e_L2_down = 0;
for ee = 1 : n_el
    for qua = 1 : n_int
        u_h   = 0;
        x_l   = 0;
        x_ele = x_coor(IEN(ee,:));
        for aa = 1 : n_en
            u_h = u_h + disp(IEN(ee,aa)) * PolyShape(pp,aa,xi(qua),0);
            x_l = x_l + x_ele(aa) * PolyShape(pp,aa,xi(qua),0);
        end
        e_L2_up   = e_L2_up + weight(qua) * (u_h - u(x_l))^2;
        e_L2_down = e_L2_down + weight(qua) * u(x_l)^2;
    end
end
e_L2_up   = sqrt(e_L2_up);
e_L2_down = sqrt(e_L2_down);
e_L2      = e_L2_up/e_L2_down;


%e_H1
e_H1_up   = 0;
e_H1_down = 0;
for ee = 1 : n_el
    for qua = 1 : n_int
        u_h_x  = 0;
        x_l    = 0;
        x_ele  = x_coor(IEN(ee,:));
        dx_dxi = 0.0;

        for aa = 1 : n_en
            x_l = x_l + x_ele(aa) * PolyShape(pp,aa,xi(qua),0);
            dx_dxi = dx_dxi + x_ele(aa) * PolyShape(pp, aa, xi(qua), 1);
        end
        dxi_dx = 1.0 / dx_dxi;

        for aa = 1 : n_en
            u_h_x = u_h_x + disp(IEN(ee,aa)) * PolyShape(pp,aa,xi(qua),1) * dxi_dx;
        end
        
        e_H1_up   = e_H1_up + weight(qua) * (u_h_x - u_x(x_l))^2;
        e_H1_down = e_H1_down + weight(qua) * u_x(x_l)^2;
    end
end
e_H1_up   = sqrt(e_H1_up);
e_H1_down = sqrt(e_H1_down);
e_H1      = e_H1_up/e_H1_down;



end

