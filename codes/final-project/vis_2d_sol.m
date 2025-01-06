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
clear all;
clc;
clf;

% 假设加载的 displacement 是一个 25 x 2 的矩阵
load("U.mat");  % 加载 displacement 数据，假设文件中包含 displacement 变量

% 网格尺寸
n_el_x = 4;  % 4 个单元 (5 个节点)
n_el_y = 4;  % 4 个单元 (5 个节点)

% 计算网格步长
hh_x = 1.0 / n_el_x;
hh_y = 1.0 / n_el_y;

% 节点数目
n_np_x = n_el_x + 1;
n_np_y = n_el_y + 1;

% 生成网格坐标
[X, Y] = meshgrid(0 : hh_x : 1, 0 : hh_y : 1);

% 从 displacement 中提取位移数据
u_x = displacement(:, 1)*1e9;  % x 方向的位移
u_y = displacement(:, 2)*1e9;  % y 方向的位移

% 将位移应用到原始网格节点
X_new = X + reshape(u_x, n_np_x, n_np_y)';
Y_new = Y + reshape(u_y, n_np_x, n_np_y)';

% 绘制变形后的网格
figure;
hold on;

% 绘制变形后的节点
plot(X_new, Y_new, 'ko', 'MarkerFaceColor', 'r');  % 红色的变形后节点

% 绘制变形后的网格线
for i = 1:n_np_y
    plot(X_new(i, :), Y_new(i, :), 'b');  % 绘制每一行
end
for j = 1:n_np_x
    plot(X_new(:, j), Y_new(:, j), 'b');  % 绘制每一列
end

% 原始网格（如果需要显示）
plot(X, Y, 'k--');  % 使用虚线绘制原始网格

% 设置图形
axis equal;
xlabel('x');
ylabel('y');
title('变形后的网格');
hold off;

% EOF