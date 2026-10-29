function P = pair_alignment(rate, side, pair, word, wordRate, nNull)
%PAIR_ALIGNMENT Leave-one-pair-out cosine alignment of word-pair difference vectors, against random word pairs.
%   rate: occurrences x units; side, pair, word: analogy labels of each occurrence; wordRate: mean rate of every
%   word type of the stories (word types x units), from which the random word pairs are drawn
if nargin < 6, nNull = 10000; end
rng(0, 'twister');
rows = [];
U = zeros(0, size(rate, 2));
for r = unique(pair(pair > 0))'
    a = slot(rate, word, pair == r & side == -1);
    b = slot(rate, word, pair == r & side == 1);
    if isempty(a) || isempty(b), continue; end
    d = b / norm(b) - a / norm(a);
    U(end + 1, :) = d / norm(d); %#ok<AGROW>
    rows(end + 1, 1) = r; %#ok<AGROW>
end
n = numel(rows);

T = NaN(size(wordRate));
for k = 1:size(wordRate, 1)
    if norm(wordRate(k, :)) > 0, T(k, :) = wordRate(k, :) / norm(wordRate(k, :)); end
end
ok = find(all(isfinite(T), 2));
R = NaN(nNull, n, size(rate, 2));
for s = 1:nNull
    m = 0;
    while m < n
        ia = ok(randi(numel(ok)));
        ib = ok(randi(numel(ok)));
        if ia == ib, continue; end
        d = T(ib, :) - T(ia, :);
        if norm(d) > 0
            m = m + 1;
            R(s, m, :) = d / norm(d);
        end
    end
end

sim = zeros(n, 1); p = zeros(n, 1); null = zeros(nNull, n);
Rfirst = permute(R(:, 1:n-1, :), [2 3 1]);
for r = 1:n
    sim(r) = mean(U([1:r-1, r+1:n], :) * U(r, :)');
    null(:, r) = squeeze(mean(pagemtimes(Rfirst, U(r, :)'), 1));
    p(r) = (1 + sum(null(:, r) >= sim(r))) / (1 + nNull);
end
P = struct('rows', rows, 'diff', U, 'sim', sim, 'null', null, 'p', p, 'q', bh_fdr(p));
end

function v = slot(rate, word, in)
w = unique(word(in));
V = zeros(numel(w), size(rate, 2));
for k = 1:numel(w)
    V(k, :) = mean(rate(in & word == w(k), :), 1);
end
v = (ones(1, numel(w)) / numel(w)) * V;
end
