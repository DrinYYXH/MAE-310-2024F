function [e0,e1] = Error(mesh,displacement,n_int,weight,xi,eta,u)

%unpack mesh
IEN = mesh.IEN;
ID  = mesh.ID;

IDH        = mesh.IDH;
ID_abandon = mesh.ID_abandon;

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

%unpack u
exact    = u.exact;
exact_x  = u.exact_x;
exact_y  = u.exact_y;
exact_xx = u.exact_xx;
exact_yy = u.exact_yy;
exact_xy = u.exact_xy;
exact_yx = u.exact_yx;



Error_Sn0 = 0;
Error_Sn1 = 0;
Exact_Sn2 = 0;


for ee = 1 : n_el
    x_ele   = x_coor( IEN(ee, 1:n_en) );
    y_ele   = y_coor( IEN(ee, 1:n_en) );
    error_Sn0 = [0,0];
    error_Sn1 = [0,0];
    exact_Sn2 = [0,0];

    for ll = 1 : n_int

            x_l     = 0.0;
            y_l     = 0.0;
            dx_dxi  = 0.0; dx_deta = 0.0;
            dy_dxi  = 0.0; dy_deta = 0.0;


            for aa = 1 : n_en
                x_l = x_l + x_ele(aa) * Quad(aa, xi(ll), eta(ll));
                y_l = y_l + y_ele(aa) * Quad(aa, xi(ll), eta(ll));
    
                [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
                dx_dxi  = dx_dxi  + x_ele(aa) * Na_xi;
                dx_deta = dx_deta + x_ele(aa) * Na_eta;
                dy_dxi  = dy_dxi  + y_ele(aa) * Na_xi;
                dy_deta = dy_deta + y_ele(aa) * Na_eta;
            end

            detJ = dx_dxi * dy_deta - dx_deta * dy_dxi;
        
         for i = 1 : n_sd  
             u_h     = 0;
             u_h_x   = 0;
             u_h_y   = 0;
            for aa = 1 : n_en
    
                [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
                Na_x = (Na_xi * dy_deta - Na_eta * dy_dxi) / detJ;
                Na_y = (-Na_xi * dx_deta + Na_eta * dx_dxi) / detJ;
          
                u_h   = u_h    + displacement(IEN(ee,aa),i) * Quad(aa, xi(ll), eta(ll));
                u_h_x = u_h_x  + displacement(IEN(ee,aa),i) * Na_x;
                u_h_y = u_h_y  + displacement(IEN(ee,aa),i) * Na_y;

            end
                error_Sn0(i) = error_Sn0(i) + weight(ll) * (u_h - exact(x_l , y_l, i))^2 * detJ;
    
                error_Sn1(i) = error_Sn1(i) + weight(ll) * ( ( u_h   - exact(x_l , y_l, i)   )^2 ...
                                                           + ( u_h_x - exact_x(x_l , y_l, i) )^2 ...
                                                           + ( u_h_y - exact_y(x_l , y_l, i) )^2 ) * detJ;
    
                exact_Sn2(i) = exact_Sn2(i) + weight(ll) * (  exact(x_l , y_l, i)    ^2 ...
                                                           +  exact_x(x_l , y_l, i)  ^2 ...
                                                           +  exact_y(x_l , y_l, i)  ^2 ...
                                                           +  exact_xx(x_l , y_l, i) ^2 ...
                                                           +  exact_yy(x_l , y_l, i) ^2 ...
                                                           +  exact_xy(x_l , y_l, i) ^2 ...
                                                           +  exact_yx(x_l , y_l, i) ^2) * detJ;

            
        end
    end
    disp(error_Sn0);
    disp(error_Sn1);
    disp(exact_Sn2);

    Error_Sn0 = Error_Sn0 + sum(error_Sn0);
    Error_Sn1 = Error_Sn1 + sum(error_Sn1);
    Exact_Sn2 = Exact_Sn2 + sum(exact_Sn2);

end
  

Error_Sn0 = sqrt(Error_Sn0);
Error_Sn1 = sqrt(Error_Sn1);
Exact_Sn2 = sqrt(Exact_Sn2);

e0 = Error_Sn0/Exact_Sn2;
e1 = Error_Sn1/Exact_Sn2;



end

