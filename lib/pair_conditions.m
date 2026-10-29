function cond = pair_conditions(side, pair, rows)
%PAIR_CONDITIONS Occurrences of each side of each word pair, ordered A1, B1, A2, B2, ...
cond = false(numel(side), 2 * numel(rows));
for r = 1:numel(rows)
    cond(:, 2 * r - 1) = pair == rows(r) & side == -1;
    cond(:, 2 * r) = pair == rows(r) & side == 1;
end
end
