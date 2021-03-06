%% system parameters

% The switched affine dynamics are defined by
%       x+ = A{sigma} * x + B{sigma} * u + E{sigma} * w + f{sigma}
N = 8; % number of steps
d = 3; % number of UAVs
n = 2 * d; % dimension of state space
m = 1; % dimension of input space
l = d; % dimension of disturbance space
ns = 2; % number of switching modes

Kh = 0.25;
Kv = 1;

topo = cell(ns,1);
topo{1} = [0, 1, 0;
           0, 0, 1;
           1, 0, 0];
topo{2} = [0, 0, 1;
           1, 0, 0;
           0, 1, 0];
%topo{2} = topo{1};

%{
%orig
Mh = 10;
offset = 0;
Mv = 5;
uH = 2;
wH = 0.1;
%}

%new
Mh = 25;
offset = 0;
Mv = 10;
uH = 5;
wH = 0.1;

A = cell(ns,1); % A matrices for each mode
for mode = 1:ns
    degrees = (sum(topo{mode}, 1));
    Aol = kron(eye(d), [1,1;0,1]);
    Acl = zeros(size(Aol));
    for i = 1:d
        scale = 1 / (degrees(i) + 1);
        for j = 1:d
            if j == i
                Acl(2*i, 2*i-1) = scale * -Kh * degrees(i);
                Acl(2*i, 2*i)   = scale * -Kv * degrees(i);
            else
                Acl(2*i, 2*j-1) =  scale * topo{mode}(j,i) * Kh;
                Acl(2*i, 2*j)   =  scale * topo{mode}(j,i) * Kv;
            end
        end
    end
    
    A{mode} = Aol + Acl;
end

B = cell(ns,1); % B matrices for each mode
for i = 1:ns
    B{i} = [0;1;zeros(n-2, 1)];
    %B{i} = repmat([0;1], d,1);
    %B{i} = kron(eye(d), [0;1]);
end

E = cell(ns,1); % E matrices for each mode
for i = 1:ns
    E{i} = kron(eye(d),[0;1]);
    %E{i} = zeros(n,0);
end

f = cell(ns,1); % f vectors for each mode
for i = 1:ns
    f{i} = zeros(n,1);
end


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

Aconv = zeros(d * d, n);
for i = 1:d
    for j = 1:d
        if i ~= j
            Aconv((i-1) * d + j, 2*i) = 1;
            Aconv((i-1) * d + j, 2*j) = -1;
        end
    end
end
bconv = ones(size(Aconv, 1), 1);
near_polyv = Polyhedron(Aconv, bconv);

Apos = [-kron(eye(d), [1, 0]);
        kron(eye(d), [0, 1]);
        -kron(eye(d), [0, 1])];
bpos = [Mh * ones(d, 1);
        Mv * ones(d,1);
        Mv * ones(d,1)];

%X = (kron(eye(d), [Mh, 0; 0, Mv]) * Polyhedron.unitBox(n)) &  (delta * near_poly); % the safe state space
%X = kron(eye(d), [Mh, 0; 0, Mv]) * Polyhedron.unitBox(n);
X = Polyhedron(Apos, bpos);% & kron(eye(d), [Mh, 0; 0, Mv]) * Polyhedron.unitBox(n);
U = uH * Polyhedron.unitBox(m); % the input space
W = wH * Polyhedron.unitBox(l); % the disturbance space

Omega =  offset*repmat([1;0],d,1) + (kron(eye(d), [1, 0; 0, 0.25]) * Polyhedron.unitBox(n)) & (0.25 * near_polyv) & (1 * near_poly); % the seed set
%Omega = (0.04 * near_polyv) & (0.12 * near_poly); % the seed set

pslsys = PolySwitchLinSys(A,X,B,U,E,W,f); % class representing the system

yalmipOptions = sdpsettings('verbose', 1); % options for the LP solver

%% j

% 4 8
multset = [4];
satind = zeros(size(multset));
elap_time = zeros(size(multset));
for i = 1:length(multset)
    start_time = cputime;
    mult = multset(i);
    % determine if the seed set generates an invariant set with respect to the
    % system dynamics
    [diagnostics, sequences, Kx_map, Kw_map, uc_map] = evalSwitchedInvariance(mult * Omega, pslsys, N, yalmipOptions);
    if diagnostics.problem == 0
        satind(i) = 1;
    end
    elap_time(i) = cputime - start_time;
end

[min(multset(find(satind))), max(multset(find(satind)))]

if diagnostics.problem == 0
    %disp('Omega generates an invariant set');
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


%% struct

st.N = N; % number of steps
st.d = d; % number of UAVs
st.n = n; % dimension of state space
st.m = m; % dimension of input space
st.l = l; % dimension of disturbance space
st.ns = ns; % number of switching modes

st.Kh = Kh;
st.Kv = Kv;

st.topo = topo;

st.uH = uH;
st.wH = wH;

st.pslsys = pslsys;

st.Omega = Omega;
st.scale = scale;
st.multset = multset;
st.satind = satind;
st.diagnostics = diagnostics;
st.elap_time = elap_time;

%% save

filename = "data/consensus/switch/1.mat";
%save(filename, 'st');