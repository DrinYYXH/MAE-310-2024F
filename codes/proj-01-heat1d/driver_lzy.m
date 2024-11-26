% clean the memory and screen
clear; 
clc; 


% Define the external source or force and boundary data
f = @(x) -20*x.^3; % f(x) = x
g = 1.0;           % u    = g  at x = 1
h = 0.0;           % -u,x = h  at x = 0
u = @(x) x.^5;      % test solve
u_x = @(x) 5*x.^4;

% Setup the mesh
pp   = 3;              % polynomial degree
n_en = pp + 1;         % number of element or local nodes
n_el = 4;              % number of elements
n_np = n_el * pp + 1;  % number of nodal points
n_eq = n_np - 1;       % number of equations
n_int = 6;


%% mesh

%initial the mesh for equally spaced nodes
mesh = GenerateMesh(n_el,n_en);

IEN    = mesh.IEN;
ID     = mesh.ID;
x_coor = mesh.coor;
hh     = mesh.hh;


%% FEM

%solve FEM problem
displacement = FEM(mesh,h,g,n_int,f);

% filename = 'displacement.txt';
% fileID = fopen(filename, 'a');
% fprintf(fileID, '%-6.4e ', displacement);
% fprintf(fileID, '\n');
% fclose(fileID);



%% visualization

% maybe a stupid way that make a large array just call it IP array, about the same as IEN array
% IP array is mainly the find every point of uu from partial coor to total coor

tt = 4;                                 %divide the points in each partial element
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

        uu(IP(ee,xip)) = uu(IP(ee,xip)) + displacement(IEN(ee,aa)) * PolyShape(pp,aa,xii,0);

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
plot(x_xii,u(x_xii),'k');



%% error

n_int = 10; % can make the Guass digree different with FEM
[e_L2,e_H1] = Error(mesh,displacement,n_int,u,u_x);





%% h-error visualization

n_ele = [2,4,6,8,10,12,14,16];
e_L2 = zeros(1,length(n_ele));
e_H1 = zeros(1,length(n_ele));


for i = 1 : length(n_ele)

    % Define the external source or force and boundary data
    f = @(x) -20*x.^3; % f(x) = x
    g = 1.0;           % u    = g  at x = 1
    h = 0.0;           % -u,x = h  at x = 0
    u = @(x) x.^5;      % test solve
    u_x = @(x) 5*x.^4;

    % Setup the mesh
    pp   = 3;              % polynomial degree
    n_en = pp + 1;         % number of element or local nodes
    n_el = n_ele(i);              % number of elements
    n_np = n_el * pp + 1;  % number of nodal points
    n_eq = n_np - 1;       % number of equations
    n_int = 10;

    %initial the mesh for equally spaced nodes
    mesh = GenerateMesh(n_el,n_en);

    IEN    = mesh.IEN;
    ID     = mesh.ID;
    x_coor = mesh.coor;
    hh     = mesh.hh;

    %solve FEM problem
    displacement = FEM(mesh,h,g,n_int,f);

    n_int = 20;
    [e_L2(i),e_H1(i)] = Error(mesh,displacement,n_int,u,u_x);


end




figure(2)
plot(log(1./n_ele),log(e_L2),'-xr');
hold on;
plot(log(1./n_ele),log(e_H1),'-.k');

xlabel('log(h)');
ylabel('log(error)');
title('h-error');
legend('e_{L2}', 'e_{H1}','Location','best');



p1 = polyfit(log(1./n_ele),log(e_L2),1);
p2 = polyfit(log(1./n_ele),log(e_H1),1);


slope1 = sprintf('Slope e_{H1} = %.2f', p2(1));
slope2 = sprintf('Slope e_{L2} = %.2f', p1(1));

% 根据需要调整文字位置
text(-2, -5.5, slope1, 'FontSize', 12, 'Color', 'k');
text(-2, -7, slope2, 'FontSize', 12, 'Color', 'r');










% EOF