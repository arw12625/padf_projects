%% system parameters

% The switched affine dynamics are defined by
%       x+ = A{sigma} * x + B{sigma} * u + E{sigma} * w + f{sigma}
N = 5; % number of steps
n = 5; % dimension of state space
m = 1; % dimension of input space
l = n; % dimension of disturbance space
ns = 2; % number of switching modes

topo = cell(ns,1);
topo{1} = [0, 1, 0, 0, 1;
           1, 0, 1, 0, 0;
           0, 1, 0, 1, 0;
           0, 0, 1, 0, 1;
           0, 0, 0, 1, 0];
topo{2} = [0, 1, 1, 1, 1;
           1, 0, 0, 0, 0;
           1, 0, 0, 0, 0;
           1, 0, 0, 0, 0;
           1, 0, 0, 0, 0];

A = cell(ns,1); % A matrices for each mode
for i = 1:ns
    A{i} = diag(1 ./ (sum(topo{i}, 2) + 1)) * (eye(n) + topo{i});
end

B = cell(ns,1); % B matrices for each mode
for i = 1:ns
    B{i} = ones(n, 1);
end

E = cell(ns,1); % E matrices for each mode
for i = 1:ns
    E{i} = eye(n);
end

f = cell(ns,1); % f vectors for each mode
for i = 1:ns
    f{i} = zeros(n,1);
end

Xm = 10;
Um = 10;
Wm = 0.1825;
delta = 2;

Acon = zeros(n * (n-1), n);
for i = 1:n
    for j = 1:n
        Acon((i-1) * n + j, i) = 1;
        Acon((i-1) * n + j, j) = -1;
    end
end
bcon = delta * ones(size(Acon, 1), 1);

near_poly = Polyhedron(Acon, bcon);

X = (Xm * Polyhedron.unitBox(n)) & near_poly; % the safe state space
U = Um * Polyhedron.unitBox(m); % the input space
W = (Wm * Polyhedron.unitBox(l)); % the disturbance space

Omega = (0.5 * Xm * Polyhedron.unitBox(n)) & (0.75 * near_poly); % the seed set

pslsys = PolySwitchLinSys(A,X,B,U,E,W,f); % class representing the system

yalmipOptions = sdpsettings('verbose', 1); % options for the LP solver

%%

% determine if the seed set generates an invariant set with respect to the
% system dynamics
[diagnostics, sequences, Kx_map, Kw_map, uc_map] = evalSwitchedInvariance(Omega, pslsys, N, yalmipOptions);


if diagnostics.problem == 0
    disp('Omega generates an invariant set');
    %{
    diagnostics.problem
    for i = 1:(size(sequences, 1) - 1)
        len_seq = sequences{i};
        for j = 1:size(len_seq, 1)
            seq = len_seq{j}
            value(Kx_map(seq))
            value(Kw_map(seq))
        end
    end
    %}
end