%% Matlab mesh
%% input, Created by Gmsh
%% ASCII
clear msh;
msh.nbNod = 9;
msh.POS = [
0 0 0;
1 0 0;
1 1 0;
0 1 0;
0.499999999998694 0 0;
1 0.499999999998694 0;
0.5000000000020591 1 0;
0 0.5000000000020591 0;
0.5000000000003766 0.5000000000003766 0;
];
msh.MAX = max(msh.POS);
msh.MIN = min(msh.POS);
msh.LINES =[
 1 5 1
 5 2 1
 2 6 2
 6 3 2
 3 7 3
 7 4 3
 4 8 4
 8 1 4
];
msh.QUADS =[
 1 5 9 8 5
 8 9 7 4 5
 5 2 6 9 5
 9 6 3 7 5
];
