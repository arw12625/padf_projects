%% system parameters

% The switched affine dynamics are defined by
%       x+ = A{sigma} * x + B{sigma} * u + E{sigma} * w + f{sigma}
N = 20; % number of steps
d = 2 * 2; % number of UAVs
n = 2 * d; % dimension of state space
m = 1; % dimension of input space
l = d; % dimension of disturbance space
ns = 1; % number of switching modes

Kh = 0.25;
Kv = 1;

topo = createGridGraph(sqrt(d), sqrt(d));

%{
g = digraph(topo{1});
bins = conncomp(g, 'Type', 'weak');
isConnected = all(bins == 1)
%}

Mh = 250;
offset = 0;
Mv = 100;
uH = 25;
%Mh = 25;
%offset = 0;
%Mv = 10;
%uH = 5;
wH = 0.1;

degrees = sum(topo, 1);
Aol = kron(eye(d), [1,1;0,1]);
Acl = zeros(size(Aol));
for i = 1:d
    scale = 1 / (degrees(i) + 1);
    for j = 1:d
        if j == i
            Acl(2*i, 2*i-1) = scale * -Kh * degrees(i);
            Acl(2*i, 2*i)   = scale * -Kv * degrees(i);
        else
            Acl(2*i, 2*j-1) =  scale * topo(j,i) * Kh;
            Acl(2*i, 2*j)   =  scale * topo(j,i) * Kv;
        end
    end
end

A = Aol + Acl;

%B{i} = [0;1;zeros(n-2, 1)];
B = [zeros(floor(n/2)-1, 1); 0;1;zeros(ceil(n/2)-1, 1)];
%B{i} = repmat([0;1], d,1);
%B{i} = kron(eye(d), [0;1]);


E = kron(eye(d),[0;1]);
%E{i} = zeros(n,0);

f = zeros(n,1);



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

sys = LTVSSys.constructLTISys(N,A,B,E,f,X * U,X,W);

yalmipOptions = sdpsettings('verbose', 1); % options for the LP solver


%% j

multset = 1;
satind = zeros(size(multset));
elap_time = zeros(size(multset));
for i = 1:length(multset)
    start_time = cputime;
    mult = multset(i);
    % determine if the seed set generates an invariant set with respect to the
    % system dynamics
    [is_rec, controller] = testAffineRecurrence(sys,Omega,N);
    if is_rec
        satind(i) = 1;
    end
    elap_time(i) = cputime - start_time;
end

[min(multset(find(satind))), max(multset(find(satind)))]



