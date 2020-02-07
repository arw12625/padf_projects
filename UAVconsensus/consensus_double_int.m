%% system parameters

% The switched affine dynamics are defined by
%       x+ = A{sigma} * x + B{sigma} * u + E{sigma} * w + f{sigma}
N = 6; % number of steps
d = 5; % number of UAVs
n = 2 * d; % dimension of state space
m = 1; % dimension of input space
l = d; % dimension of disturbance space
ns = 2; % number of switching modes

Kh = 0.5;
Kv = 0.9;

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
for mode = 1:ns
    degrees = (sum(topo{mode}, 1));
    Aol = kron(eye(d), [1,1;0,0.99]);
    Acl = zeros(size(Aol));
    for i = 2:d
        scale = 1 / (degrees(i) + 1);
        for j = 1:d
            if j == i
                Acl(2*i, 2*i-1) = scale * -Kh * degrees(i);
                Acl(2*i, 2*i)   = scale * -Kv * degrees(i);
            else
                Acl(2*i, 2*j-1) =  scale * topo{mode}(i,j) * Kh;
                Acl(2*i, 2*j)   =  scale * topo{mode}(i,j) * Kv;
            end
        end
    end
    
    A{mode} = Aol + Acl;
end

B = cell(ns,1); % B matrices for each mode
for i = 1:ns
    B{i} = [0;1;zeros(n-2, 1)];
end

E = cell(ns,1); % E matrices for each mode
for i = 1:ns
    E{i} = kron(eye(d),[0;1]);
end

f = cell(ns,1); % f vectors for each mode
for i = 1:ns
    f{i} = zeros(n,1);
end

Mh = 10;
Mv = 10;
delta = 10;
uH = 10;
wH = 0.01;

Acon = zeros(d * d, n);
for i = 1:d
    for j = 1:d
        if i ~= j
            Acon((i-1) * d + j, 2*i-1) = 1;
            Acon((i-1) * d + j, 2*j-1) = -1;
        end
    end
end
bcon = ones(size(Acon, 1), 1);
near_poly = Polyhedron(Acon, bcon);

X = (kron(eye(d), [Mh, 0; 0, Mv]) * Polyhedron.unitBox(n)) &  (delta * near_poly); % the safe state space
U = uH * Polyhedron.unitBox(m); % the input space
W = wH * Polyhedron.unitBox(l); % the disturbance space

Omega = (kron(eye(d), [0.5, 0; 0, 0.05]) * Polyhedron.unitBox(n)) & (0.375*near_poly); % the seed set

pslsys = PolySwitchLinSys(A,X,B,U,E,W,f); % class representing the system

yalmipOptions = sdpsettings('verbose', 1); % options for the LP solver

%% j

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