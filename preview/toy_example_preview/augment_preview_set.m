function [preview_Omega] = augment_preview_set(plsys, Omega, preview)

if preview == 0
    preview_Omega = Omega;
else

    preview_A = blkdiag(Omega.A, kron(eye(preview), plsys.W.A));
    preview_b = [Omega.b; repmat(plsys.W.b, preview, 1)];

    preview_Omega = Polyhedron(preview_A, preview_b);
end
end

