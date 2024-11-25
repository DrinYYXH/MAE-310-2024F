function [ID,IEN] = GenerateMesh(n_el,n_en)
%MESH 此处显示有关此函数的摘要
%   此处显示详细说明
n_np = n_el * (n_en - 1) + 1;

IEN = zeros(n_el, n_en);

for ee = 1 : n_el
  for aa = 1 : n_en
    IEN(ee, aa) = (ee - 1) * (n_en - 1) + aa;
  end
end

% Setup the ID array for the problem
ID = 1 : n_np;
ID(end) = 0;% this depend on the B.C.


end

