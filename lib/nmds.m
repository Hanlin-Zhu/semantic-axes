function [Y, stress, D] = nmds(M)
%NMDS Nonmetric 3-D MDS of the rows of M (conditions x units), after z-scoring each unit, on cosine distances.
Z = zscore(M);
D = pdist2(Z, Z, 'cosine');
D = (D + D') / 2;
D = D - diag(diag(D));
[Y, stress] = mdscale(D, 3, 'Criterion', 'stress');
end
