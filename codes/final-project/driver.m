clear;
clc;
 
%% User-defined Parts
D_defined = 1;  % define the asumption for two-dimentional analysis plane stain(1) or plane stress(2)








%% auto-going

% define Const.
E      = 1e9;
nu     = 0.3;
lambda = E*nu / ( (1 + nu) * (1 - 2*nu) );
miu    = E / (2*(1 + nu));

% stiffness matrix
if D_defined == 1
    D = [lambda + 2*miu , lambda         , 0
         lambda         , lambda + 2*miu , 0
         0              , 0              , miu];
elseif D_defined == 2
    D = [1              , nu             , 0
         nu             , 1              , 0
         0              , 0              , (1 - nu)/2] * (E / (1 - nu^2));
else
    error('Error: value of D_defined should be 1 or 2');
end






















%EOF