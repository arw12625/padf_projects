%% System definitions

% System matrices
% x - state
% u - input
% w - disturbance
% x(t+1) = A x(t) + B u(t) + E w(t) + f
n = 8;
m = 4;
l = n;
wscale = 0.05;

% The horizon considered
T = 8;

A = normrnd(0,1,n,n);
B = normrnd(0,1,n,m);
E = eye(n);%normrnd(0,1,n,l);
f = zeros(n,1);

% System constraints as polyhedra
Xsafe = Polyhedron.unitBox(n);
Usafe = Polyhedron.unitBox(m);
XU = Xsafe * Usafe;
Xterm = Xsafe;
W = wscale * Polyhedron.unitBox(l);

% Construct a linear time-varying switched system (LTVSS) representing the linear
% time-invariant double integrator system.
% For now all methods operate on LTVSS to avoid redundant functionality
sys = LTVSSys.constructLTISys(T,A,B,E,f,XU,Xterm,W);

%%

% The seed set used for computation
Omegabase = Polyhedron.unitBox(n);

scaleRange = [0.2];

tic
for scale = scaleRange
% Determine if there is an affine disturbance feedback controller that can
% drive the system starting in Omega into Omega in exactly T steps subject
% to all disturbances and constraints.
% If the system is recurrent (satisfies this condition), then the method
% returns such a controller.
    Omega = scale * Omegabase;
    Omega = outerInvar;
    % create a copy of the system to add target constraint
    csys = copy(sys);
    horizon = T;
    
    
    % bisect?
    is_rec = zeros(size(Omega.H,1),1);
    for i = 1:size(Omega.H, 1)
        
        Omega_face = Polyhedron(...
            'H', Omega.H, ...
            'He', [Omega.He; Omega.H(i,:)]); 
        csys.Xterm = Omega;
        
        initialConditions = computeInitialConditions(csys, Omega_face, horizon);
        admissibleTrajectories = computeLiftedAdmissibleTrajectories(csys, horizon);

        yalmipOptions = sdpsettings('verbose', 1, 'solver', ''); % options for the LP solver
        [diagnostics, affineController] = computeAffineController(csys, initialConditions, admissibleTrajectories, yalmipOptions);
        
        is_rec(i) = (diagnostics.problem == 0);
        if ~is_rec(i)
            break
        end
    end
    
    
    %[is_rec, controller] = testAffineRecurrence(sys,Omega,T);
    if all(is_rec)
        found = 1;
        break
    end
end
toc
rank(ctrb(A,B)) == n


%%
tic
preUnion = computePreUnion(sys, Omega, T);
toc

%%

% Compute an outer approximation of the maximal invariant set by
% iteratively applying pre to the safe set.
tic
outerInvar = computeOuterInvar(sys, Xsafe, T);
toc
