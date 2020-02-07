is_rec_inv = 1
while is_rec_inv == 1
%% System definitions

found = 0;
while found == 0
% System matrices
% x - state
% u - input
% w - disturbance
% x(t+1) = A x(t) + B u(t) + E w(t) + f
n = 2;
A = rand(n);
B = rand(n);
E = eye(n);
f = zeros(n,1);

% System constraints as polyhedra
Xsafe = Polyhedron.unitBox(n);
Usafe = Polyhedron.unitBox(n);
XU = Xsafe * Usafe;
Xterm = Xsafe;
wscale = rand(1);
W = wscale * Polyhedron.unitBox(n);

% The horizon considered
T = 10;

% Construct a linear time-varying switched system (LTVSS) representing the linear
% time-invariant double integrator system.
% For now all methods operate on LTVSS to avoid redundant functionality
sys = LTVSSys.constructLTISys(T,A,B,E,f,XU,Xterm,W);

% The seed set used for computation
Omegabase = Polyhedron.unitBox(n);

for scale = 0.1:0.1:1
% Determine if there is an affine disturbance feedback controller that can
% drive the system starting in Omega into Omega in exactly T steps subject
% to all disturbances and constraints.
% If the system is recurrent (satisfies this condition), then the method
% returns such a controller.
Omega = scale * Omegabase;
[is_rec, controller] = testAffineRecurrence(sys,Omega,T);
if is_rec == 1
    found = 1;
    break
end
end
end
%%

% Compute an outer approximation of the maximal invariant set by
% iteratively applying pre to the safe set.
%tic
outerInvar = computeOuterInvar(sys, Xsafe, 8);
%toc

%%

insys = LTVSSys.constructLTISys(1,A,B,E,f,XU,outerInvar,W);
[is_rec_inv, controller_inv] = testAffineRecurrence(insys,outerInvar,1);

end