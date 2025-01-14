1. 如何运行计算过程
    进行计算过程的之程序为driver.m，它将主要调用网格生成器MeshGenerate.m，计算函数FEM.m进行计算。
driver.m包括上下两个运行区块，上半部分为单独计算某一特定网格尺寸的程序，下半部风为计算不同网格尺寸
并绘制误差分析图的程序。

2. 如何设置网格尺寸
    在单一计算中，网格尺寸可以通过直接修改driver.m中38、39行内容来设定网格尺寸。这样设定的尺寸会通
过结构体自动传输至其他运算单元。
    这一操作仅限于使用自定义网格生成器有效，如果是调用gmesh生成网格，可自动忽略上述两个位置的输入。

3. 如何调用自定义网格生成器
    自定义网格生成器位于MeshGenerate.m的下半区块，请注释掉上半区块内容，并保证下半区块内容未被注释。

4. 如何调用gmesh进行网格划分
    gmesh网格生成器位于MeshGenerate.m的上半区块，请注释掉下半区块内容，并保证上半区块内容未被注释。
    请注意以下几点：
        1）确保gmesh网格生成器中调用gmesh_output.m的步骤被注释，并且在脚本文件gmesh_input.m中设定好
           需求，之后主动运行该脚本。
        2）确保gmesh网格生成器中调用gmesh_output.m的步骤被注释，因为该脚本还有bug未修复。
        3）谁当并运行gmesh_output.m后需要使用gmsh.exe打开生成的脚本文件input.geo，直接进行导出，设定
           文件名为output.m，的MATLAB版本，并取消勾选“Save all elements”。

5.如何修改边界条件
    边界条件的修改位置在MeshGenerate.m，在指定的网格生成器模块下找到Dirichlet_BC矩阵，其为一个4*2的矩
    阵，第一行的每一个数字风别表示了x方向的边界条件指针，第二行为y方向。1、2、3、4分别对应下、右、上、左
    四条边。
    若希望设定某一条边的某一个方向为Dirichlet条件，则需要将其对应数字位置的数字生定为其对应数字，反若希
望是Neumann边界条件则设为0。

6.如何可视化
    可视化程序为vis_2d_sol.m，在运行了计算程序后直接运行该程序即可。

7.注意事项
    1）修改边界条件需要在driver.m中手动输入边界条件的方程（f或是h）。
    2）源代码中提供了两个范例与driver.m，请选择一个运行，另外一个保持注释状态。
    3）整个程序仅能够使用四边形网格。