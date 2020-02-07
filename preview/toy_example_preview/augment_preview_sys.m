function [preview_pslsys] = augment_preview_sys(pslsys, preview)

if preview == 0
    preview_pslsys = pslsys;
else
    shift_mat = [eye(preview-1); zeros(1,preview-1)];

    Ap = cell(1);
    Bp = cell(1);
    Ep = cell(1);
    fp = cell(1);
    for i = 1:pslsys.ns
        Ap{i} = blkdiag([pslsys.A{i}, pslsys.E{i}], kron(shift_mat, eye(pslsys.l))); 
        Bp{i} = [pslsys.B{i}; zeros(pslsys.l * preview, pslsys.m)];
        Ep{i} = [zeros(pslsys.n + pslsys.l * (preview-1), pslsys.l); eye(pslsys.l)];
        fp{i} = [pslsys.f{i}; zeros(pslsys.l * preview, 1)];
    end
    nCon = size(pslsys.X.A,1);
    Xp = Polyhedron([pslsys.X.A, zeros(nCon, pslsys.l * preview)], pslsys.X.b);
    Up = pslsys.U;
    Wp = pslsys.W;

    preview_pslsys = PolySwitchLinSys(Ap, Xp, Bp, Up, Ep, Wp, fp);
end
end

