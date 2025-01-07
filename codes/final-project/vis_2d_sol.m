clear all; clc; clf;

load("U.mat");

hh_x = 1.0 / n_el_x;
hh_y = 1.0 / n_el_y;

n_np_x = n_el_x + 1;
n_np_y = n_el_y + 1;

[X, Y] = meshgrid(0 : hh_x : 1, 0 : hh_y : 1);
Z = reshape(displacement, n_np_x, n_np_y)';
surf(X, Y, Z);

shading interp

az = -61;
el = 20;
view(az,el);

% EOF


%%
clf;

% 假设加载的 displacement 是一个 25 x 2 的矩阵
load("U.mat");  % 加载 displacement 数据，假设文件中包含 displacement 变量

% 网格尺寸
n_el_x = 3;  % 4 个单元 (5 个节点)
n_el_y = 3;  % 4 个单元 (5 个节点)

% 计算网格步长
hh_x = 1.0 / n_el_x;
hh_y = 1.0 / n_el_y;

% 节点数目
n_np_x = n_el_x + 1;
n_np_y = n_el_y + 1;

% 生成网格坐标
[X, Y] = meshgrid(0 : hh_x : 1, 0 : hh_y : 1);

% 从 displacement 中提取位移数据
u_x_num = displacement(:, 1)*1e9;  % x 方向的位移
u_y_num = displacement(:, 2)*1e9;  % y 方向的位移

% 将数值解位移应用到原始网格节点
X_num = X + reshape(u_x_num, n_np_x, n_np_y)';
Y_num = Y + reshape(u_y_num, n_np_x, n_np_y)';

% 绘制数值解变形后的网格
figure(1);
hold on;

% 绘制数值解变形后的节点
plot(X_num, Y_num, 'ko', 'MarkerFaceColor', 'r');  % 红色的变形后节点

% 绘制数值解变形后的网格线
for i = 1:n_np_y
    plot(X_num(i, :), Y_num(i, :), 'b');  % 绘制每一行
end
for j = 1:n_np_x
    plot(X_num(:, j), Y_num(:, j), 'b');  % 绘制每一列
end

% 原始网格（绘制四条边框）
for i = 1:n_np_y
    plot(X(i, :), Y(i, :), 'k--');  % 绘制每一行的虚线
end
for j = 1:n_np_x
    plot(X(:, j), Y(:, j), 'k--');  % 绘制每一列的虚线
end

% 设置图形
axis equal;
xlabel('x');
ylabel('y');
title('数值解变形后的网格');
hold off;



% 计算理论解的位移（新增部分）
u_x_ext = zeros(n_np_x, n_np_y);  % 初始化理论解 x 方向位移矩阵
u_y_ext = zeros(n_np_x, n_np_y);  % 初始化理论解 y 方向位移矩阵

% 从 exact 中提取位移数据
for i = 1:n_np_x
    for j = 1:n_np_y
        % 计算每个网格点 (x, y) 上的理论位移
        u_x_ext(i,j) = exact(X(i,j), Y(i,j), 1)*1e9;  % x 方向的理论位移
        u_y_ext(i,j) = exact(X(i,j), Y(i,j), 2)*1e9;  % y 方向的理论位移
    end
end

% 将理论解位移应用到原始网格节点
X_ext = X + reshape(u_x_ext, n_np_x, n_np_y)';
Y_ext = Y + reshape(u_y_ext, n_np_x, n_np_y)';

% 绘制理论解变形后的网格
figure(2);
hold on;

% 绘制理论解变形后的节点
plot(X_ext, Y_ext, 'ko', 'MarkerFaceColor', 'r');  % 红色的变形后节点

% 绘制理论解变形后的网格线
for i = 1:n_np_y
    plot(X_ext(i, :), Y_ext(i, :), 'b');  % 绘制每一行
end
for j = 1:n_np_x
    plot(X_ext(:, j), Y_ext(:, j), 'b');  % 绘制每一列
end

% 原始网格（绘制四条边框）
for i = 1:n_np_y
    plot(X(i, :), Y(i, :), 'k--');  % 绘制每一行的虚线
end
for j = 1:n_np_x
    plot(X(:, j), Y(:, j), 'k--');  % 绘制每一列的虚线
end

% 设置图形
axis equal;
xlabel('x');
ylabel('y');
title('理论解变形后的网格');
hold off;


% EOF