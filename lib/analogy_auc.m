function A = analogy_auc(rate, side, word, nPerm)
%ANALOGY_AUC Trial-level and word-level balanced AUC of every unit for one analogy, with permutation P values.
%   rate: occurrences x units; side: -1 or +1 for occurrences in the analogy, 0 otherwise; word: word ID
if nargin < 4, nPerm = 1000; end
rng(0, 'twister');
use = side ~= 0;
X = rate(use, :);
y = double(side(use));
[~, ~, wi] = unique(word(use));
wordSign = zeros(max(wi), 1);
wordSign(wi) = y;
W = zeros(max(wi), size(X, 2));
for j = 1:size(X, 2)
    W(:, j) = accumarray(wi, X(:, j)) ./ accumarray(wi, 1);
end

nU = size(rate, 2);
A = struct('trial', nan(nU, 1), 'word', nan(nU, 1), 'p_trial', nan(nU, 1), 'p_word', nan(nU, 1), ...
    'p_within', nan(nU, 1));
idxA = find(y == -1);
idxB = find(y == 1);
tokA = find(wordSign == -1);
tokB = find(wordSign == 1);
[~, posA] = ismember(wi(idxA), tokA);
[~, posB] = ismember(wi(idxB), tokB);
nA = accumarray(posA, 1, [numel(tokA) 1]);
nB = accumarray(posB, 1, [numel(tokB) 1]);
assignA = repelem((1:numel(tokA))', nA);
assignB = repelem((1:numel(tokB))', nB);
nPos = sum(y == 1);
nPosW = sum(wordSign == 1);
isB = [false(numel(tokA), 1); true(numel(tokB), 1)];
for j = 1:nU
    r = tiedrank(X(:, j));
    A.trial(j) = balanced(r, y == 1);
    rw = tiedrank(W(:, j));
    A.word(j) = balanced(rw, wordSign == 1);
    A.p_trial(j) = perm_p(r, nPos, A.trial(j), nPerm);
    A.p_word(j) = perm_p(rw, nPosW, A.word(j), nPerm);
    xA = X(idxA, j);
    xB = X(idxB, j);
    null = zeros(nPerm, 1);
    for b = 1:nPerm
        rateA = accumarray(assignA, xA(randperm(numel(idxA)))) ./ nA;
        rateB = accumarray(assignB, xB(randperm(numel(idxB)))) ./ nB;
        null(b) = balanced(tiedrank([rateA; rateB]), isB);
    end
    A.p_within(j) = (1 + sum(null >= A.word(j))) / (nPerm + 1);
end
end

function a = balanced(ranks, pos)
n1 = sum(pos);
n0 = numel(pos) - n1;
a = (sum(ranks(pos)) - n1 * (n1 + 1) / 2) / (n1 * n0);
a = max(a, 1 - a);
end

function p = perm_p(ranks, nPos, obs, nPerm)
n = numel(ranks);
null = zeros(nPerm, 1);
for b = 1:nPerm
    q = randperm(n);
    null(b) = (sum(ranks(q(1:nPos))) - nPos * (nPos + 1) / 2) / (nPos * (n - nPos));
end
null = max(null, 1 - null);
p = (1 + sum(null >= obs)) / (nPerm + 1);
end
