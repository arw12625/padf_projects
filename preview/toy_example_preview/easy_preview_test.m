
preview = 1;

A = {[1,1,0;0,0,1;0,0,0]};
B = {[0;0;1]};
E = {[1;0;0]};
f = {[0;0;0]};

X = 8 * Polyhedron.unitBox(3);
U = 2 * Polyhedron.unitBox(1);
W = 1 * Polyhedron.unitBox(1);

pslsys = PolySwitchLinSys(A,X,B,U,E,W,f);
prev_sys = augment_preview_sys(pslsys, preview);

Omega = 1 * Polyhedron.unitBox(3);
prev_Omega = augment_preview_set(pslsys, Omega, preview);

N = 6;

options = sdpsettings('verbose', 1);

[diagnostics, sequences, Kx_map, Kw_map, uc_map] = evalSwitchedInvariance(prev_Omega, prev_sys, N, options);

gen_invar = 1 - diagnostics.problem;

save("preview/toy_example_preview/models/red", "pslsys", "preview", "N", "gen_invar")