function acc = ccgp(conds, nIter, null, nWorkers)
%CCGP Cross-condition generalization performance, leave-one-pair-out, with an L2 logistic decoder.
%   conds: 1 x 2n cell of trials x units firing rates, ordered A1, B1, A2, B2, ...
%   null: 'none', 'neuron' (unit identities shuffled per condition) or 'label' (class labels shuffled)
%   acc: nIter x 1 mean accuracy over held-out pairs
if nargin < 4, nWorkers = 0; end
nU = size(conds{1}, 2);
nTr = min(cellfun(@(x) size(x, 1), conds));
pos = 2:2:numel(conds);
neg = 1:2:numel(conds);
grid = [1 / nU 3e-3 1e-2 3e-2 6e-2 1e-1 3e-1];
acc = NaN(nIter, 1);
parfor (it = 1:nIter, nWorkers)
    rng(it, 'threefry');
    c = conds;
    if strcmp(null, 'neuron')
        c = shuffle_units(c);
    end
    perf = NaN(numel(pos), 1);
    for k = 1:numel(pos)
        trP = pos([1:k-1, k+1:end]);
        trN = neg([1:k-1, k+1:end]);
        if numel(trP) < 2
            lam = 1 / nU;
        else
            lam = pick_lambda(c, trP, trN, nTr, grid);
        end
        [Xtr, Xte, ytr, yte] = sample_trials(c, trP, trN, pos(k), neg(k), nTr);
        if strcmp(null, 'label')
            ytr = ytr(randperm(numel(ytr)));
            yte = yte(randperm(numel(yte)));
        end
        Xtr = zscore(Xtr);
        Xte = zscore(Xte);
        keep = ~any(isnan(Xtr), 1) & ~any(isnan(Xte), 1);
        if ~any(keep), continue; end
        mdl = fitclinear(Xtr(:, keep), ytr, 'Lambda', lam, 'Learner', 'logistic', 'Regularization', 'ridge');
        perf(k) = mean(predict(mdl, Xte(:, keep)) == yte);
    end
    acc(it) = mean(perf, 'omitnan');
end
end

function lam = pick_lambda(c, trP, trN, nTr, grid)
err = NaN(numel(grid), 1);
for g = 1:numel(grid)
    e = NaN(numel(trP), 1);
    for j = 1:numel(trP)
        [Xtr, Xva, ytr, yva] = sample_trials(c, setdiff(trP, trP(j)), setdiff(trN, trN(j)), trP(j), trN(j), nTr);
        Xtr = zscore(Xtr);
        Xva = zscore(Xva);
        keep = ~any(isnan(Xtr), 1) & ~any(isnan(Xva), 1);
        if ~any(keep), continue; end
        mdl = fitclinear(Xtr(:, keep), ytr, 'Lambda', grid(g), 'Learner', 'logistic', 'Regularization', 'ridge');
        e(j) = mean(predict(mdl, Xva(:, keep)) ~= yva);
    end
    err(g) = mean(e, 'omitnan');
end
if all(isnan(err))
    lam = exp(mean(log(grid)));
    return
end
best = find(abs(err - min(err)) < 1e-12);
lam = exp(mean(log(grid(best))));
end

function [Xtr, Xte, ytr, yte] = sample_trials(c, trP, trN, teP, teN, nTr)
nU = size(c{1}, 2);
Xtr = zeros(nTr * (numel(trP) + numel(trN)), nU);
Xte = zeros(2 * nTr, nU);
for u = 1:nU
    Xtr(:, u) = [draw(c, trP, u, nTr); draw(c, trN, u, nTr)];
    Xte(:, u) = [draw(c, teP, u, nTr); draw(c, teN, u, nTr)];
end
ytr = [ones(nTr * numel(trP), 1); -ones(nTr * numel(trN), 1)];
yte = [ones(nTr, 1); -ones(nTr, 1)];
end

function v = draw(c, set, u, nTr)
v = zeros(nTr * numel(set), 1);
for s = 1:numel(set)
    x = c{set(s)}(:, u);
    v((s - 1) * nTr + (1:nTr)) = x(randperm(numel(x), nTr));
end
end

function c = shuffle_units(c)
for i = 1:numel(c)
    M = c{i}(:, 1);
    for u = 2:size(c{i}, 2)
        n = min(size(M, 1), size(c{i}, 1));
        M = M(randperm(size(M, 1), n), :);
        x = c{i}(randperm(size(c{i}, 1), n), u);
        M(:, u) = x;
    end
    c{i} = M;
end
for i = 1:numel(c)
    c{i} = c{i}(:, randperm(size(c{i}, 2)));
end
end
