% 文件路径
geo_file = 'input.geo';

% 打开文件，准备写入
fid = fopen(geo_file, 'w');

% 编写几何内容
fprintf(fid, '// 定义几何体：一个矩形区域\n');
fprintf(fid, '\n');
fprintf(fid, '// dh = 1.0\n');
fprintf(fid, '\n');
fprintf(fid, 'Point(1) = {0, 0, 0, dh}; // 点1\n');
fprintf(fid, 'Point(2) = {1, 0, 0, dh}; // 点2\n');
fprintf(fid, 'Point(3) = {1, 1, 0, dh}; // 点3\n');
fprintf(fid, 'Point(4) = {0, 1, 0, dh}; // 点4\n');
fprintf(fid, '\n');
fprintf(fid, 'Line(1) = {1, 2}; // 连接点1和点2的线段\n');
fprintf(fid, 'Line(2) = {2, 3}; // 连接点2和点3的线段\n');
fprintf(fid, 'Line(3) = {3, 4}; // 连接点3和点4的线段\n');
fprintf(fid, 'Line(4) = {4, 1}; // 连接点4和点1的线段\n');
fprintf(fid, '\n');
fprintf(fid, 'Line Loop(1) = {1, 2, 3, 4}; // 将四条线段构成一个线圈\n');
fprintf(fid, 'Plane Surface(1) = {1}; // 创建面1\n');
fprintf(fid, '\n');
fprintf(fid, 'Mesh 2; // 生成二维网格\n');
fprintf(fid, '\n');
fprintf(fid, 'Physical Line("Bottom") = {1}; // 命名Bottom边\n');
fprintf(fid, 'Physical Line("Right") = {2}; // 命名Right边\n');
fprintf(fid, 'Physical Line("Top") = {3}; // 命名Top边\n');
fprintf(fid, 'Physical Line("Left") = {4}; // 命名Left边\n');
fprintf(fid, '\n');
fprintf(fid, 'Physical Surface("Domain") = {1}; // 设置整个面为物理区域\n');
fprintf(fid, 'Recombine Surface{1}; //将三角形网格改为四边形网格');

% 关闭文件
fclose(fid);

disp('Geo文件已生成');


