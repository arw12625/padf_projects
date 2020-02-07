%% Base System

A = {[1,1,0;0,0,1;0,0,0]};
B = {[0;0;1]};
E = {[1;0;0]};
f = {[0;0;0]};

X = 8 * Polyhedron.unitBox(3);
U = 2 * Polyhedron.unitBox(1);
W = 1 * Polyhedron.unitBox(1);
pslsys = PolySwitchLinSys(A,X,B,U,E,W,f);

%% Preview = 1

prev_sys1 = augment_preview_sys(pslsys, 1);

pinv1 = computeOuterApproxInvariantSwitch(prev_sys1, 8);

%% Preview = 2

Omega = 1 * Polyhedron.unitBox(3);
preview = 2;
prev_sys2 = augment_preview_sys(pslsys, preview);
prev_Omega = augment_preview_set(pslsys, Omega, preview);

N = 6;

options = sdpsettings('verbose', 1);

[diagnostics, sequences, Kx_map, Kw_map, uc_map] = evalSwitchedInvariance(prev_Omega, prev_sys2, N, options);
gen_invar = 1 - diagnostics.problem

pinv2 = computeConvInvar(prev_sys2, prev_Omega, N);

%maximal
%pinv2 = computeOuterApproxInvariantSwitch(prev_sys2, 8);

%% Compare

% For preview = 1, we compute the slice corresponding to w1=0
a1 = slice(pinv1,[4]);

% For preview = 2, we compute the slice corresponding to w1=0 and then the
% projection removing w2
a2 = projection(slice(pinv2,[4]),[1,2,3]);

a1 < a2
volume(a1)
volume(a2)

plot(a1)
plot(a2)

%%

prod1 = pinv1 * W;
prod1 < pinv2
volume(prod1)
volume(pinv2)
