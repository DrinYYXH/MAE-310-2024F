clf;

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
n_sd   = mesh.n_sd;

n_ee   = n_en * n_sd;

x_coor = mesh.x_coor;
y_coor = mesh.y_coor;

hx = mesh.hx;
hy = mesh.hy;


load("U.mat");  % 加载 displacement 数据，假设文件中包含 displacement 变量


% 从 displacement 中提取位移数据
u_x_num = displacement(:, 1);  % x 方向的位移
u_y_num = displacement(:, 2);  % y 方向的位移

% 将数值解位移应用到原始网格节点
X_num = x_coor + u_x_num;
Y_num = y_coor + u_y_num;

% 绘制数值解变形后的网格
figure(1);
hold on;

% 绘制数值解变形后的节点
plot(X_num, Y_num, 'ko', 'MarkerFaceColor', 'r','MarkerSize',2);  % 红色的变形后节点

% % 绘制数值解变形后的网格线
% plot(X_num, Y_num, 'b');

% % 原始网格（绘制四条边框）
% for ny = 1 : n_np_x
%     % 绘制竖直方向上的虚线
%     plot([x_coor(ny), x_coor(ny)], [y_coor(1), y_coor(n_np_y)], 'k--');
% end
% for nx = 1 : n_np_y
%     % 绘制水平方向上的虚线
%     plot([x_coor(1), x_coor(n_np_x)], [y_coor(nx), y_coor(nx)], 'k--');
% end

% 设置图形
axis equal;
xlabel('x');
ylabel('y');
title('数值解变形后的网格');
hold off;




% 计算理论解的位移（新增部分）
u_x_ext = zeros(n_np,1);  % 初始化理论解 x 方向位移矩阵
u_y_ext = zeros(n_np,1);  % 初始化理论解 y 方向位移矩阵

% 从 exact 中提取位移数据
for ny = 1 : n_np_y
    for nx = 1 : n_np_x
        index = (ny-1)*n_np_x + nx;
        % 计算每个网格点 (x, y) 上的理论位移
        u_x_ext(index) = exact(x_coor(index), y_coor(index), 1);  % x 方向的理论位移
        u_y_ext(index) = exact(x_coor(index), y_coor(index), 2);  % y 方向的理论位移
    end
end
% 将理论解位移应用到原始网格节点
x_ext = x_coor + u_x_ext;
y_ext = y_coor + u_y_ext;

% 绘制理论解变形后的网格
figure(2);
hold on;

% 绘制理论解变形后的节点
plot(x_ext, y_ext, 'ko', 'MarkerFaceColor', 'r','MarkerSize',2);  % 红色的变形后节点

% 绘制理论解变形后的网格线
for ny = 1:n_np_y
    for nx = 1:n_np_x
    plot([x_ext(1), x_ext(n_np_x)], [y_ext(ny), y_ext(ny)], 'b');  % 绘制每一行
    plot([x_ext(nx), x_ext(nx)], [y_ext(1), y_ext(n_np_y)], 'b');  % 绘制每一列
    end
end

% % 原始网格（绘制四条边框）
% for ny = 1 : n_np_x
%     % 绘制竖直方向上的虚线
%     plot([x_coor(ny), x_coor(ny)], [y_coor(1), y_coor(n_np_y)], 'k--');
% end
% for nx = 1 : n_np_y
%     % 绘制水平方向上的虚线
%     plot([x_coor(1), x_coor(n_np_x)], [y_coor(nx), y_coor(nx)], 'k--');
% end

% 设置图形
axis equal;
xlabel('x');
ylabel('y');
title('理论解变形后的网格');
hold off;


% EOF