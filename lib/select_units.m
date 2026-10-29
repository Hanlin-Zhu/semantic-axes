function u = select_units(A, cutoff)
%SELECT_UNITS Units with trial-level p(AUC) below cutoff (at least the 8 lowest), plus units with both word-level
%   P values below 0.01.
u = find(A.p_trial < cutoff);
if numel(u) < 8
    [~, u] = mink(A.p_trial, 8);
end
u = union(u, find(A.p_within < .01 & A.p_word < .01));
end
