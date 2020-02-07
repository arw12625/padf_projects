clc;clear all; close all;

param.preview = 2;

param.u_max = 2;
param.u_min = -2;
param.d_max = 2;
param.d_min = -2;

dyn_aug = get_1d_dyn_aug(param);

Safe = Polyhedron('lb', -32, 'ub', 32);

P = dyn_aug.D;

Safe_aug = Safe;

for i = 1:param.preview
    Safe_aug = Safe_aug*P;
end


[W_aug,log] = dyn_aug.win_always(Safe_aug,0,1,1);
