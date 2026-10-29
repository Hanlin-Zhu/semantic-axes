function S = fig2(g, auc, nIter, nWorkers, out, azel)
%FIG2 Gender word pairs in the population: MDS (A), tail-aligned directions (B), alignment of each pair with the
%   other pairs (C) and cross-condition generalization (D), with the tuned units of the gender analogy.
units = select_units(auc, 0.125);
X = g.rate(:, units);
P = pair_alignment(X, g.side, g.pair, g.word, g.word_rate(:, units));
cond = pair_conditions(g.side, g.pair, P.rows);
n = numel(P.rows);
M = zeros(2 * n, numel(units));
for c = 1:2 * n
    M(c, :) = mean(X(cond(:, c), :), 1);
end
[Y, stress] = nmds(M);
sig = P.q < .05;
C = analogy_colors();
C([1 8], :) = C([8 1], :);
C = C(1:n, :);
la = g.labels(P.rows, 1);
lb = g.labels(P.rows, 2);

f = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 16 13]);
plot_pairs(gca, Y, la, lb, C, sig);
view(azel);
title(sprintf('Gender (%d units); Kruskal stress = %.3f', numel(units), stress), 'FontWeight', 'normal');
exportgraphics(f, fullfile(out, 'Fig2A.pdf'), 'ContentType', 'image', 'Resolution', 300);

f = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 10 10]);
plot_globe(gca, Y(2:2:end, :) - Y(1:2:end, :), C, sig);
exportgraphics(f, fullfile(out, 'Fig2B.pdf'), 'ContentType', 'image', 'Resolution', 300);

f = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 9 n * 0.8 + 2]);
plot_alignment(gca, P, C, lb + " - " + la);
title('Gender', 'FontWeight', 'normal', 'HorizontalAlignment', 'right');
exportgraphics(f, fullfile(out, 'Fig2C.pdf'), 'ContentType', 'vector');

conds = arrayfun(@(c) X(cond(:, c), :), 1:2 * n, 'UniformOutput', false);
acc = ccgp(conds, nIter, 'none', nWorkers);
nullNeuron = ccgp(conds, nIter, 'neuron', nWorkers);
nullLabel = ccgp(conds, nIter, 'label', nWorkers);
obs = mean(acc);
pLabel = (1 + sum(nullLabel >= obs)) / (nIter + 1);
pNeuron = (1 + sum(nullNeuron >= obs)) / (nIter + 1);
f = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 6 5]);
plot_ccgp(gca, obs, {nullLabel, nullNeuron}, [pLabel pNeuron]);
exportgraphics(f, fullfile(out, 'Fig2D.pdf'), 'ContentType', 'vector');

writetable(table(lb + " - " + la, P.sim, P.p, P.q, 'VariableNames', {'pair', 'cosine', 'p', 'q'}), ...
    fullfile(out, 'Fig2C.csv'));
S = struct('units', units, 'alignment', P, 'mds', Y, 'stress', stress, 'ccgp', acc, 'null_label', nullLabel, ...
    'null_neuron', nullNeuron, 'ccgp_mean', obs, 'p_label', pLabel, 'p_neuron', pNeuron);
save(fullfile(out, 'Fig2.mat'), '-struct', 'S');
fprintf('Fig 2: %d units; %d/%d pairs aligned; Kruskal stress %.3f; CCGP %.3f (label P = %.4f, neuron P = %.4f)\n', ...
    numel(units), nnz(sig), n, stress, obs, pLabel, pNeuron);
end

function plot_pairs(ax, Y, la, lb, C, sig)
hold(ax, 'on');
grey = [0.8 0.8 0.8];
for r = [find(~sig); find(sig)]'
    a = Y(2 * r - 1, :);
    b = Y(2 * r, :);
    c = C(r, :);
    if sig(r)
        plot3(ax, [a(1) b(1)], [a(2) b(2)], [a(3) b(3)], 'Color', [c 0.8], 'LineWidth', 3);
        arrow(ax, a, b, [0.35 0.35 0.35], '-', 1);
    else
        c = grey;
        arrow(ax, a, b, grey, '--', 1);
    end
    scatter3(ax, [a(1) b(1)], [a(2) b(2)], [a(3) b(3)], 75, c, 'filled', 'MarkerEdgeColor', [0.2 0.2 0.2], ...
        'MarkerFaceAlpha', 0.9, 'MarkerEdgeAlpha', 0.6);
    text(ax, a(1), a(2), a(3), " " + la(r), 'FontName', 'Courier', 'FontSize', 14, 'Color', c, 'Interpreter', 'none');
    text(ax, b(1), b(2), b(3), " " + lb(r), 'FontName', 'Courier', 'FontSize', 14, 'Color', c, 'Interpreter', 'none');
end
A = Y(1:2:end, :);
B = Y(2:2:end, :);
u = B(sig, :) - A(sig, :);
u = mean(u ./ vecnorm(u, 2, 2), 1);
mL = mean(A(sig, :), 1);
mR = mL + u / norm(u) * norm(mean(B(sig, :), 1) - mL);
arrow(ax, mL - 0.1 * (mR - mL), mR + 0.1 * (mR - mL), [0 0 1], '-', 3);
text(ax, mL(1), mL(2), mL(3), ' mean direction', 'Color', 'b', 'FontSize', 14);
axis(ax, 'tight');
xlabel(ax, 'MDS-1'); ylabel(ax, 'MDS-2'); zlabel(ax, 'MDS-3');
set(ax, 'FontSize', 12, 'TickLength', [0.01 0.01]);
xl = xlim(ax); yl = ylim(ax); zl = zlim(ax);
h = patch(ax, xl([1 2 2 1]), yl([1 1 2 2]), zl(1) * [1 1 1 1], [0.9 0.9 0.9], 'FaceAlpha', 0.5, 'EdgeColor', 'none');
uistack(h, 'bottom');
end

function arrow(ax, a, b, c, style, width)
quiver3(ax, a(1), a(2), a(3), b(1) - a(1), b(2) - a(2), b(3) - a(3), 0, 'Color', c, 'LineStyle', style, ...
    'LineWidth', width, 'MaxHeadSize', min(0.6, 0.12 / norm(b - a)));
end

function plot_globe(ax, V, C, sig)
V = V ./ vecnorm(V, 2, 2);
mu = mean(V(sig, :), 1);
mu = mu / norm(mu);
h = [0 0 1];
if abs(dot(h, mu)) > 0.95, h = [0 1 0]; end
u0 = h - dot(h, mu) * mu;
u0 = u0 / norm(u0);
v0 = cross(mu, u0);
Q = [V * u0', V * v0'];
Q = Q(sig, :) - mean(Q(sig, :), 1);
[E, L] = eig(Q' * Q / max(size(Q, 1) - 1, 1));
[~, k] = max(diag(L));
pc = E(1, k) * u0 + E(2, k) * v0;
R = rotation(mu, [0 0 1]);
V = V * R';
pc = (R * (pc / norm(pc))')';
r = 0.8;
[sx, sy, sz] = sphere(80);
surf(ax, r * sx, r * sy, r * sz, 'FaceAlpha', 0.04, 'EdgeAlpha', 0.06, 'FaceColor', [0.67 0.67 0.67]);
hold(ax, 'on');
C(~sig, :) = repmat([0.7 0.7 0.7], nnz(~sig), 1);
for i = 1:size(V, 1)
    quiver3(ax, 0, 0, 0, r * V(i, 1), r * V(i, 2), r * V(i, 3), 0, 'Color', C(i, :), 'LineWidth', 1.6, 'MaxHeadSize', 0.25);
    scatter3(ax, r * V(i, 1), r * V(i, 2), r * V(i, 3), 16, C(i, :), 'filled', 'MarkerFaceAlpha', 0.7, 'MarkerEdgeColor', 'none');
end
quiver3(ax, 0, 0, 0, 0, 0, 1, 0, 'Color', 'b', 'LineWidth', 1.5, 'MaxHeadSize', 0.15);
text(ax, 0, 0, 1.12, 'mean direction', 'Color', 'b', 'FontSize', 12, 'HorizontalAlignment', 'center');
axis(ax, 'equal', 'off');
xlim(ax, [-1.05 1.05]); ylim(ax, [-1.05 1.05]); zlim(ax, [-0.8 1.15]);
camtarget(ax, [0 0 0]); camup(ax, [0 0 1]); campos(ax, 3.5 * pc);
end

function R = rotation(a, b)
v = cross(a, b);
c = dot(a, b);
s = norm(v);
if s < 1e-12
    R = eye(3) * sign(c);
    return
end
K = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
R = eye(3) + K + K * K * ((1 - c) / s ^ 2);
end

function plot_alignment(ax, P, C, names)
n = numel(P.sim);
s = 0.8;
sig = P.q < .05;
hold(ax, 'on');
for r = 1:n
    [f, x] = ksdensity(P.null(:, r), 'NumPoints', 200);
    f = f / max(f) * 0.38;
    fill(ax, [x fliplr(x)], [s * r - f, s * r * ones(size(f))], 0.85 * [1 1 1], 'EdgeColor', 'none', 'FaceAlpha', 0.7);
end
C(~sig, :) = 0.5;
scatter(ax, P.sim, s * (1:n), 36, C, 'filled');
plot(ax, [0 0], [0 n], 'k');
lims = [min([P.null(:); P.sim]) max([P.null(:); P.sim])];
for r = find(sig)'
    text(ax, P.sim(r) + 0.05 * diff(lims), s * r + 0.15, '*', 'FontSize', 12, 'FontWeight', 'bold', ...
        'VerticalAlignment', 'middle');
end
names(~sig) = "\color[rgb]{0.7 0.7 0.7}" + names(~sig);
set(ax, 'YTick', s * (1:n), 'YTickLabel', names, 'YDir', 'reverse', 'TickDir', 'out', 'TickLength', [0.003 0.003], ...
    'FontSize', 12, 'Box', 'off', 'Layer', 'top');
ax.YAxis.FontName = 'Courier';
ylim(ax, s * [0.4 n + 0.4]);
xlim(ax, [-0.17, lims(2) + 0.06 * diff(lims)]);
xticks(ax, [0 0.3]);
xlabel(ax, 'Cosine similarity');
end

function plot_ccgp(ax, obs, nulls, p)
pos = 0.75 * [1 2];
hold(ax, 'on');
for i = 1:2
    [f, x] = ksdensity(nulls{i}, 'NumPoints', 200);
    f = f / max(f) * 0.35;
    fill(ax, [pos(i) - f, pos(i) - f(1)], [x, x(1)], 0.85 * [1 1 1], 'EdgeColor', 'none', 'FaceAlpha', 0.7);
end
plot(ax, pos, [obs obs], 'r');
scatter(ax, pos, [obs obs], 36, 'r', 's', 'filled');
plot(ax, [0.3 1.95], [0.5 0.5], 'k--');
for i = find(p < .05)
    text(ax, pos(i), obs + 0.02, '*', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'FontSize', 12, 'FontWeight', 'bold');
end
set(ax, 'XTick', [0.65 1.45], 'XTickLabel', {'label\newlinenull', 'neuron\newlinenull'}, 'YTick', [0 0.5 1], ...
    'TickDir', 'out', 'Box', 'off', 'FontSize', 13);
xlim(ax, [0.3 1.75]); ylim(ax, [0 1.15]); ylabel(ax, 'CCGP');
end
