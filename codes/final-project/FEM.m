function displacement = FEM(mesh,n_int,weight,xi,eta,f,kappa)

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

% allocate the stiffness matrix and load vector
K = spalloc(n_eq, n_eq, 9 * n_eq);
F = zeros(n_eq, 1);


% loop over element to assembly the matrix and vector
for ee = 1 : n_el
    x_ele = x_coor( IEN(ee, 1:n_en) );
    y_ele = y_coor( IEN(ee, 1:n_en) );

    k_ele = zeros(n_en, n_en); % element stiffness matrix
    f_ele = zeros(n_en, 1);    % element load vector
    for i = 1 : n_sd
        for ll = 1 : n_int
            x_l = 0.0; y_l = 0.0;
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
                Na = Tri(aa, xi(ll), eta(ll));
                [Na_xi, Na_eta] = Tri_grad(aa, xi(ll), eta(ll));
                Na_x = (Na_xi * dy_deta - Na_eta * dy_dxi) / detJ;
                Na_y = (-Na_xi * dx_deta + Na_eta * dx_dxi) / detJ;

                B_a = [Na_x , 0
                       0    , Na_y
                       Na_y , Na_x];

                f_ele(n_sd * (aa - 1) + i) = f_ele(n_sd * (aa - 1) + i) + weight(ll) * detJ * f(x_l, y_l, i) * Na;

                for bb = 1 : n_en
                    Nb = Tri(bb, xi(ll), eta(ll));
                    [Nb_xi, Nb_eta] = Tri_grad(bb, xi(ll), eta(ll));
                    Nb_x = (Nb_xi * dy_deta - Nb_eta * dy_dxi) / detJ;
                    Nb_y = (-Nb_xi * dx_deta + Nb_eta * dx_dxi) / detJ;

                    B_b = [Nb_x , 0
                           0    , Nb_y
                           Nb_y , Nb_x];

                    k_ele(aa, bb) = k_ele(aa,bb) + weight(ll) * detJ * B_a' * D * B_b;
                end % end of bb loop
            end % end of aa loop
        end % end of quadrature loop

        for aa = 1 : n_en
            PP = ID(i,IEN(ee,aa));
            if PP > 0
                F(PP) = F(PP) + f_ele(n_sd * (aa - 1) + i);

                for bb = 1 : n_en
                    QQ = ID(i,IEN(ee,bb));
                    if QQ > 0
                        K(PP, QQ) = K(PP, QQ) + k_ele(aa, bb);
                    else
                        % modify F with the boundary data
                        % here we do nothing because the boundary data g is zero or
                        % homogeneous
                    end
                end
            end
        end
    end
end


% solve the stiffness matrix
dn = K \ F;

% insert dn back into the vector for all nodes
displacement = zeros(n_np, 1);

for ii = 1 : n_np
  index = ID(ii);
  if index > 0
    displacement(ii) = dn(index);
  else
    % modify disp with the g data. Here it does nothing because g is zero
  end
end

% save the solution vector and number of elements to disp with name
% HEAT.mat
save("HEAT", "displacement", "n_el_x", "n_el_y");
end

