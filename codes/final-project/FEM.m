function displacement = FEM(mesh,n_int,weight,xi,eta,f,D,g)

%unpack mesh
IEN = mesh.IEN;
ID  = mesh.ID;
LM  = mesh.LM;

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

% allocate the stiffness matrix and load vector
K = spalloc(n_eq, n_eq, 9 * n_eq);
F = zeros(n_eq, 1);

e{1} = [1; 0];
e{2} = [0; 1];

% loop over element to assembly the matrix and vector
for ee = 1 : n_el
    x_ele = x_coor( IEN(ee, 1:n_en) );
    y_ele = y_coor( IEN(ee, 1:n_en) );
    
    k_ele = zeros(n_ee, n_ee); % element stiffness matrix
    f_ele = zeros(n_ee, 1);    % element load vector
    f_g   = zeros(n_ee, 1);    % Dirichlet B.C.
    f_h   = zeros(n_ee, 1);    % Neumann B.C.
    
     for ll = 1 : n_int
            x_l = 0.0; y_l = 0.0;
            dx_dxi = 0.0; dx_deta = 0.0;
            dy_dxi = 0.0; dy_deta = 0.0;
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

          for aa = 1 : n_en
                Na = Quad(aa, xi(ll), eta(ll));
                [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
                Na_x = (Na_xi * dy_deta - Na_eta * dy_dxi) / detJ;
                Na_y = (-Na_xi * dx_deta + Na_eta * dx_dxi) / detJ;

                B_a = [Na_x , 0;
                       0    , Na_y;
                       Na_y , Na_x];
            for i = 1 : n_sd
                
                for bb = 1 : n_en
                    Nb = Quad(bb, xi(ll), eta(ll));
                    [Nb_xi, Nb_eta] = Quad_grad(bb, xi(ll), eta(ll));
                    Nb_x = (Nb_xi * dy_deta - Nb_eta * dy_dxi) / detJ;
                    Nb_y = (-Nb_xi * dx_deta + Nb_eta * dx_dxi) / detJ;

                    B_b = [Nb_x , 0;
                           0    , Nb_y;
                           Nb_y , Nb_x];
                    for j = 1 : n_sd
                        pp = n_sd*(aa-1) + i;
                        qq = n_sd*(bb-1) + j;
                        k_ele(pp,qq) = k_ele(pp,qq) + weight(ll) * detJ * e{i}' * B_a' * D * B_b * e{j};

                    if ID_abandon(IEN(aa,ee),i) > 0
                        f_g(qq)   = k_ele(pp,qq) * g(x_l, y_l, i);
                        f_ele(pp) = f_ele(pp) - f_g(qq);
                    end

                    end % end of j loop
                end % end of bb loop

                f_ele(pp) = f_ele(pp) + weight(ll) * detJ * f(x_l, y_l, i) * Na;

                    if IDH(IEN(aa,ee),i) > 0
                        f_ele(pp) = f_ele(pp) + Na * weight(ll) * detJ * h(x_l, y_l, i);
                    end

            end % end of i loop
        end % end of aa loop
     end


    for i = 1 : n_sd
        for aa = 1 : n_en
            pp = n_sd * (aa - 1) + i;
            PP = ID(IEN(ee,aa),i);
            if PP > 0
                F(PP) = F(PP) + f_ele(pp);
                for bb = 1 : n_en
                    qq = n_sd*(bb-1) + j;
                    QQ = ID(IEN(ee,bb),i);
                    if QQ > 0
                        for j = 1 : n_sd
                            K(PP, QQ) = K(PP, QQ) + k_ele(pp,qq);
                        end
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
displacement = zeros(n_np, n_sd);

for ii = 1 : n_np
    for jj = 1 : n_sd
        index = ID(ii,jj);
        if index > 0
            displacement(ii,jj) = dn(index);
        else
            % modify disp with the g data. Here it does nothing because g is zero
        end
    end
end
% save the solution vector and number of elements to disp with name
% HEAT.mat
save("U", "displacement", "n_el_x", "n_el_y");
end
