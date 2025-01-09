% Gmsh命令路径
gmsh_path = 'H:\CSM\gmsh-4.13.1-Windows64\gmsh-4.13.1-Windows64\gmsh.exe';  % 如果 Gmsh 已经添加到系统路径，直接使用 'gmsh'；否则，提供 Gmsh 可执行文件的绝对路径。

% 调用 Gmsh 生成 .msh 文件
command = sprintf('%s input.geo -2 -format m -o output.m', gmsh_path);
system(command);

disp('.m文件已生成');
