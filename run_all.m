% Runs the gender analysis of Figure 2A-D on the demonstration data: tuning AUCs of every unit, selection of the
% gender-tuned units, MDS, tail-aligned directions, pair alignment and CCGP. Results go to results/.

ccgpIterations = 100;
workers = 0;
mdsView = [17.0645 55.4180];

root = fileparts(mfilename('fullpath'));
addpath(fullfile(root, 'lib'));
out = fullfile(root, 'results');
if ~isfolder(out), mkdir(out); end
g = load(fullfile(root, 'data', 'gender_demo.mat'));
tic
auc = analogy_auc(g.rate, g.side, g.word);
fprintf('AUC of %d units: %.1f s\n', size(g.rate, 2), toc);
fig2(g, auc, ccgpIterations, workers, out, mdsView);
