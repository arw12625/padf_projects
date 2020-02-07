%% System definitions

% This script investigates invariant set computations for a TCL system


% Number of discrete states
K = 7;

% These adjacency matrices model the simple case where in one mode the
% discrete temperature increments by inc while in the other mode it 
% decrements by dec at each time step.
inc = 2;
dec = 1;

N0_adjacency = zeros(K,K);
N0_adjacency((1+dec):end,1:(end-dec)) = eye(K-dec);

N1_adjacency = zeros(K,K);
N1_adjacency(1:(end-inc),(inc+1):end) = eye(K-inc);

% System matrices
% x - state
% u - input
% w - disturbance
% x(t+1) = A x(t) + B u(t) + E w(t) + f
A = blkdiag(N0_adjacency, N1_adjacency);
B = [-N0_adjacency, N0_adjacency; N1_adjacency, -N1_adjacency];
E = zeros(2 * K, 0);
f = zeros(2 * K, 1);

% System constraints as polyhedra
N = 9999; % number of loads
upperBound = 3600; % upper bound on the number of loads that are on
lowerBound = 2000; % lower bound on the number of loads that are on

% This polyhedron models that constraints for the upper and lower bound on
% the number of active loads, the nonnegativity requirements, and the total
% number of loads.
powerConstraints = Polyhedron('H', ...
   [ones(1,K), zeros(1,K), N - lowerBound;
    -eye(K), zeros(K), zeros(K,1);
    zeros(1,K), ones(1,K), upperBound;
    zeros(K), -eye(K), zeros(K,1);], 'He', ...
    [ones(1,2*K), N;]);

% This polyhedron models the constraints that the number of loads to switch
% from m1 to m2 cannot be more than the number of loads at m1, and that the
% number of loads switched must be nonnegative.
jointInputConstraints = Polyhedron('H', ...
    [-eye(K), zeros(K), eye(K), zeros(K), zeros(K,1);
     zeros(K), zeros(K), -eye(K), zeros(K), zeros(K,1);
     zeros(K), -eye(K), zeros(K), eye(K), zeros(K,1);
     zeros(K), zeros(K), zeros(K), -eye(K), zeros(K,1)]);
fullInputSpace = Polyhedron('H', [zeros(1,2*K), 1]);

% This polyhedron is the combination of the state and input constraints
XU = (powerConstraints * fullInputSpace) & jointInputConstraints;
Xterm = powerConstraints;

% For now we ignore disturbance
W = Polyhedron();

% The horizon considered
T = 15;

% Construct a linear time-varying switched system (LTVSS) representing the linear
% time-invariant double integrator system.
% For now all methods operate on LTVSS to avoid redundant functionality
sys = LTVSSys.constructLTISys(T,A,B,E,f,XU,Xterm,W);

%% Invariant set computations

% The seed set used for computation
Omega = Polyhedron('H', ...
   [ones(1,K), zeros(1,K), N - lowerBound+20;
    -eye(K), zeros(K), zeros(K,1);
    zeros(1,K), ones(1,K), upperBound-20;
    zeros(K), -eye(K), zeros(K,1);
    eye(K), zeros(K), N / K * ones(K,1);
    zeros(K), eye(K), N / K * ones(K,1)], 'He', ...
    [ones(1,2*K), N;
     zeros(1,K-1),1, zeros(1,K), 0;
     zeros(1,K),1,zeros(1,K-1) 0]);
 
 
Omega = Polyhedron('H', [...
    %zeros(K-1,1), eye(K-1), zeros(K-1,K), (N-lowerBound) / (K-1)^2 * 2 * (K-1:-1:1)';
    zeros(K-1,1), eye(K-1), zeros(K-1,K), (N-lowerBound) / (K-1) * 2 * ones(K-1,1);
    %zeros(K-1,1), -eye(K-1), zeros(K-1,K), -(N-lowerBound) / (K-1) / 3 * ones(K-1,1);
    %zeros(K-1,K), eye(K-1), zeros(K-1,1), upperBound / (K-1)^2 * 2 * (1:K-1)';
    zeros(K-1,K), eye(K-1), zeros(K-1,1), upperBound / (K-1) * 2 * ones(K-1,1);
    %zeros(K-1,K), -eye(K-1), zeros(K-1,1), -upperBound / (K-1) / 3 * ones(K-1,1);
    ones(1,K), zeros(1,K), N - lowerBound * 1.3;
    zeros(1,K), ones(1,K), upperBound * 0.75], 'He', ...
    [1, zeros(1,2*K-1), 0;
     zeros(1,2*K-1), 1, 0;]);
Omega = Omega & powerConstraints;

%Omega = Polyhedron('V',[0,0,0,0,3333,3333,0,0,0,3333,0,0]);

% Determine if there is an affine disturbance feedback controller that can
% drive the system starting in Omega into Omega in exactly T steps subject
% to all disturbances and constraints.
% If the system is recurrent (satisfies this condition), then the method
% returns such a controller.
Omega.minHRep();
%[is_rec, controller] = testAffineRecurrencePiecewise(sys,Omega, Omega.getFacet(),T);
[is_rec, controller] = testAffineRecurrence(sys, Omega,T);

% Without more constraints, an invariant set is given by the entire space
% with no control.

%%

% Determine all possible reachable states using this controller over the
% horizon. If the system is recurrent, then total reach set is invariant.
[reachMap, totalReach] = reachableSetAffineController(sys, Omega, controller);

%%

% Determine the states that can reach Omega within the horizon.
% (Note technically this set computes the convex hull of the union 
% of the sets pre^i(Omega) for i = 1:T which is slightly different)
% If the system is recurrent then this set is invariant.
preUnion = computePreUnion(sys, Omega, T);

%%

% Compute an outer approximation of the maximal invariant set by
% iteratively applying pre to the safe set.
% This method has trouble with numerical stability
outerInvar = computeOuterInvar(sys, Xterm, 5);
