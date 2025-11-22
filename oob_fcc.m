function idxoob = oob_fcc(p, S, margin)
% possiblepositions: Mx3 Cartesian coords
% u1,u2,u3: 1x3 unit vectors
% S.br: half side length parameter
% S.rc: cutoff margin from boundary

% Convert all points to fractional coords
f = p(:,1:3) * S.fcc.invA;

% Inside test: all 0 <= f <= 1
inside = all(f >= 0 & f <= 1, 2);

% Distances to nearest face (in Cartesian units)
d1 = min(f(:,1), 1-f(:,1)) * norm(S.fcc.a1);
d2 = min(f(:,2), 1-f(:,2)) * norm(S.fcc.a2);
d3 = min(f(:,3), 1-f(:,3)) * norm(S.fcc.a3);
dist_to_boundary = min([d1,d2,d3],[],2);

% Margin test
margin_ok = dist_to_boundary >= margin;

% Final mask
if margin==0 
    keep_mask = inside;
elseif margin>=0
    keep_mask = inside & margin_ok;
elseif margin<0
    keep_mask=margin_ok;
end

idxoob=~keep_mask;

end