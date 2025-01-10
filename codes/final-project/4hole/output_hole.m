%% Matlab mesh
%% input, Created by Gmsh
%% ASCII
clear msh;
msh.nbNod = 16;
msh.POS = [
1 -1 0;
1 1 0;
-1 1 0;
-1 -1 0;
-0.7 -1 0;
-1 -0.7 0;
-0.7878679656440357 -0.7878679656440357 0;
-0.7228361403694678 -0.8851949699938777 0;
-0.8851949705866085 -0.7228361401239506 0;
-1 0.1499999999978956 0;
-2.750244476601438e-12 1 0;
1 2.750244476601438e-12 0;
0.1500000000021044 -1 0;
0.1060660171802601 0.1060660171802601 0;
-0.4425974852935426 0.1385819299381059 0;
0.1385819298174308 -0.4425974849944357 0;
];
msh.MAX = max(msh.POS);
msh.MIN = min(msh.POS);
msh.LINES =[
 5 8 0
 8 7 0
 7 9 0
 9 6 0
 6 10 0
 10 3 0
 3 11 0
 11 2 0
 2 12 0
 12 1 0
 1 13 0
 13 5 0
 2 14 0
 14 7 0
];
msh.QUADS =[
 3 11 15 10 0
 10 15 9 6 0
 11 2 14 15 0
 15 14 7 9 0
 2 14 16 12 0
 12 16 13 1 0
 14 7 8 16 0
 16 8 5 13 0
];
msh.PNT =[
 1 0
 2 0
 3 0
 4 0
 5 0
 6 0
 7 0
];
