R = 0.3;  //设定半径长度为0.3
L = 1.0;  //设定单位线段长度为1，应该就是整个框的边长

dx = 0.1;

Point(1) = {L, -L, 0,dx};		//设定第一个点，右下
Point(2) = {L, L, 0,dx};		//设定第二个点，右上
Point(3) = {-L, L, 0,dx};		//设定第三个点，左上
Point(4) = {-L, -L, 0,dx}; 		//设定第四个点，左下，圆心
Point(5) = {-L + R, -L, 0,dx};	//设定第五个点，右圆弧端点
Point(6) = {-L, -L + R, 0,dx};	//设定第六个点，上圆弧端点
Point(7) = {-L + Cos(Pi/4) * R, -L + Sin(Pi/4) * R, 0,dx};  //设定第七个点，圆弧中点

Circle(1) = {5, 4, 7};		//以4为圆心，从5到7画圆弧1
Circle(2) = {7, 4, 6};		//以4为圆心，从7到6画圆弧2

Line(3) = {6, 3};		//画线段3，从6到3
Line(4) = {3, 2};		//画线段4，从3到2
Line(5) = {2, 1};		//画线段5，从2到1
Line(6) = {1, 5};		//画线段6，从1到5
Line(7) = {2, 7};		//画线段7，从2到7

Curve Loop(1) = {4, 7, 2, 3};
Plane Surface(1) = {1};	//选择线段4、7、2、3，构成平面1

Curve Loop(2) = {7, -1, -6, -5};
Plane Surface(2) = {2};	//选择线段7、-1、-6、-5，构成平面1，负号代表反向走线

Transfinite Line{1, 2, 3, 4, 5, 6, 7} = 3;	//设定这几条线上各有3个网格点

Transfinite Surface{1};		//使用跨优先插值方法在目标平面构建结构化网格，使其连接上边界的节点
Transfinite Surface{2};

Recombine Surface{1};		//将三角形网格改为四边形网格
Recombine Surface{2};

Mesh.ElementOrder = 1;		//设定单元阶数为1
Mesh.Algorithm = 8;			//设定生成网格的算法

Mesh 2;

// EOF
