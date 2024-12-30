function [e0,e1] = Error(mesh,displacment,n_int,weight,xi,eta,exact,exact_x,exact_y)

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

x_coor = mesh.x_coor;
y_coor = mesh.y_coor;

hx = mesh.hx;
hy = mesh.hy;


error_Sn0 = 0;
error_Sn1 = 0;
exact_Sn1 = 0;

for ee = 1 : n_el
  x_ele   = x_coor( IEN(ee, 1:n_en) );
  y_ele   = y_coor( IEN(ee, 1:n_en) );
  u_h     = 0;
  u_h_x  = 0;
  u_h_y = 0;
  
  for ll = 1 : n_int    
    x_l = 0.0; 
    y_l = 0.0;
    dx_dxi = 0.0; dx_deta = 0.0;
    dy_dxi = 0.0; dy_deta = 0.0;
    for aa = 1 : n_en
      x_l = x_l + x_ele(aa) * Tri(aa, xi(ll), eta(ll));
      y_l = y_l + y_ele(aa) * Tri(aa, xi(ll), eta(ll));

      [Na_xi, Na_eta] = Tri_grad(aa, xi(ll), eta(ll));
      dx_dxi  = dx_dxi  + x_ele(aa) * Na_xi;
      dx_deta = dx_deta + x_ele(aa) * Na_eta;
      dy_dxi  = dy_dxi  + y_ele(aa) * Na_xi;
      dy_deta = dy_deta + y_ele(aa) * Na_eta;
    end

    detJ = dx_dxi * dy_deta - dx_deta * dy_dxi;

    for aa = 1 : n_en 

      [Na_xi, Na_eta] = Tri_grad(aa, xi(ll), eta(ll));
      Na_x = (Na_xi * dy_deta - Na_eta * dy_dxi) / detJ;
      Na_y = (-Na_xi * dx_deta + Na_eta * dx_dxi) / detJ;

      u_h   = u_h  + displacment(IEN(ee,aa)) * Tri(aa, xi(ll), eta(ll));
      u_h_x = u_h_x + displacment(IEN(ee,aa)) * Na_x;
      u_h_y = u_h_y  + displacment(IEN(ee,aa)) * Na_y;

    end


    error_Sn0 = error_Sn0 + weight(ll) * (u_h - exact(x_l , y_l))^2 * detJ;

    error_Sn1 = error_Sn1 + weight(ll) * ( ( u_h   - exact(x_l , y_l)   )^2 ...
                                         + ( u_h_x - exact_x(x_l , y_l) )^2 ...
                                         + ( u_h_y - exact_y(x_l , y_l) )^2 ) * detJ;

    exact_Sn1 = exact_Sn1 + weight(ll) * (  exact(x_l , y_l)   ^2 ...
                                         +  exact_x(x_l , y_l) ^2 ...
                                         +  exact_y(x_l , y_l) ^2 ) * detJ;
    
  end
end

error_Sn0 = sqrt(error_Sn0);
error_Sn1 = sqrt(error_Sn1);
exact_Sn1 = sqrt(exact_Sn1);

e0 = error_Sn0/exact_Sn1;
e1 = error_Sn1/exact_Sn1;



end

