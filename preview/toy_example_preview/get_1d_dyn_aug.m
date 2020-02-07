function dyn = get_1d_dyn_aug(param)

preview = param.preview;
u_max = param.u_max;
u_min = param.u_min;

d_max = param.d_max;
d_min = param.d_min;

N = 1;

constant_1d;


if preview == 0
    A_aug = [A];
else
    A_aug = [A, D, zeros(N,preview-1)
             zeros(preview-1,N+1),eye(preview-1)
             zeros(1,N+preview)];
end


B_aug = [B;zeros(preview,1)];

if preview ~= 0
    D_aug = zeros(N+preview,1);
    D_aug(end) = 1;
else
    D_aug =D;
end

K_aug = zeros(N+preview,1);

P = Polyhedron('lb', d_min, 'ub', d_max);
XU = Polyhedron('H', [zeros(1,N+preview), 1 , u_max
                      zeros(1,N+preview), -1, -u_min]);

dyn = Dyn(A_aug, K_aug, B_aug, XU, [],[],[],{zeros(N+preview)}, {D_aug}, P);    
