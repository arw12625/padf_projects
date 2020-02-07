%% System definitions

% This script demonstrates various invariant set computations for a 
% discrete time double integrator using this library. 

% System matrices
% x - state
% u - input
% w - disturbance
% x(t+1) = A x(t) + B u(t) + E w(t) + f
A = [1,1;0,0];
B = eye(2);
E = eye(2);
f = [0;0];

% System constraints as polyhedra
Xsafe = Polyhedron.unitBox(2) & Polyhedron('Ae',[1,0],'be',0);
Usafe = Polyhedron.unitBox(2) & Polyhedron('Ae',[0,1],'be',0);
XU = Xsafe * Usafe;
Xterm = Xsafe;
W = 0.25 * Polyhedron.unitBox(2) & Polyhedron('Ae',[1,0],'be',0);

% The horizon considered
T = 6;

% Construct a linear time-varying switched system (LTVSS) representing the linear
% time-invariant double integrator system.
% For now all methods operate on LTVSS to avoid redundant functionality
sys = LTVSSys.constructLTISys(T,A,B,E,f,XU,Xterm,W);

%% Invariant set computations

% The seed set used for computation
Omega = (0.4* Polyhedron.unitBox(2)) & Xsafe;

% Determine if there is an affine disturbance feedback controller that can
% drive the system starting in Omega into Omega in exactly T steps subject
% to all disturbances and constraints.
% If the system is recurrent (satisfies this condition), then the method
% returns such a controller.
[is_rec, controller] = testAffineRecurrence(sys,Omega,T);

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
%tic
outerInvar = computeOuterInvar(sys, Xsafe, 3);
%toc


%%

% Plot all of the sets computed
figure(T)
%plot([Xsafe, outerInvar, preUnion, totalReach, Omega])
%legend(["Safe Set", "Maximal Invariant", "Union of Pre", "Reach of Controller", "Seed Set"])
plot([Xsafe, outerInvar, preUnion, totalReach, Omega])
legend(["Safe Set X_s", "Maximal Invariant", "Union of Pre \Omega_N", "Reach of Controller", "Seed Set \Omega"])
title("Linear System Invariant Sets with N=6")
%{
% Plot the reachable sets
for t = 1:T
    figure(t)
    seqs = sys.sequences{1,t};
    plot(reachMap(seqs{1}))
end
%}