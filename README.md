# Semantic axes in human hippocampal neurons: gender demonstration

MATLAB code that runs the analysis behind Figure 2A–D of "A geometric foundation for word meaning in the brain" on demonstration data derived from the recordings.

The analysis covers the gender analogy:

- the tuning (AUC) of every unit;
- the selection of gender-tuned units;
- MDS of the word pairs;
- their tail-aligned directions;
- the alignment of each word pair with the other pairs;
- cross-condition generalization (CCGP).

## Demonstration data

`data/gender_demo.mat`:

- `rate` (292 × 200): the 292 occurrences of gender words in the six stories (rows, in story order) by 200 units.
  - The units are the top 200 principal components of the firing rates of the 437 recorded hippocampal neurons (91% of the variance), not single neurons.
  - Rates were counted in a window that starts 80 ms after word onset and lasts the duration of the word plus 40 ms.
- `side`: −1 or +1 for the two sides of the analogy (`sides`: male, female).
- `pair`: the gender word pair (1–12) of each occurrence.
- `word`: an anonymous ID of the word form. Singular and plural forms, for example, are different forms on the same side of a pair.
- `word_rate` (1,353 × 200): the mean of the 200 units over the occurrences of every word type of the stories, from which random word pairs are drawn. The words are not given.
- `labels` (12 × 2): pair labels, with singular and plural forms pooled. The labels are shuffled across pairs, so the label a pair carries is not its actual word pair.

The panels therefore have the form of Figure 2A–D but are not the published panels: the units are principal components, and the pair labels are shuffled. The published values are in the Source Data of the paper.

## What the code does

1. `analogy_auc`: for each unit, the balanced AUC (the AUC or 1 − AUC, whichever is larger) between the single-occurrence responses to the two sides, with P from 1,000 label shuffles. It also gives the word-level balanced AUC on each word form's mean response, with two P values: one from shuffling word labels, and one from reassigning occurrences among the word forms of each side.
2. `select_units`: the gender-tuned units. These are units with trial-level p(AUC) < 0.125 (at least the 8 lowest), plus units with both word-level P values below 0.01.
3. `pair_alignment`: for each pair, the mean cosine similarity between its normalized difference vector and those of the other pairs, against 10,000 sets of random word pairs. The P values are Benjamini–Hochberg corrected across pairs.
4. `nmds`: nonmetric 3-D MDS of the 24 word means. Each unit is z-scored across words, and distances are cosine distances. The fit is reported as Kruskal stress.
5. The tail-aligned directions are the difference vectors of the pairs in the MDS space, rotated so that the mean direction of the aligned pairs points up.
6. `ccgp`: leave-one-pair-out cross-condition generalization with an L2-regularized logistic decoder, whose λ is chosen by inner leave-one-pair-out cross-validation. It has two nulls, shuffled labels and shuffled unit identities.

## Requirements and installation

- MATLAB with the Statistics and Machine Learning Toolbox. Tested with MATLAB R2025a and Statistics and Machine Learning Toolbox 25.1 on macOS 26.7.1, on an iMac with an Apple M1 and 16 GB of memory. Any operating system that runs MATLAB should work.
- No other software and no special hardware are needed.
- The Parallel Computing Toolbox is optional; set `workers` at the top of `run_all.m` to use it.
- Installation is unzipping the folder, which takes under a minute.

## Running

In MATLAB, from this folder:

```matlab
run_all
```

`run_all` writes Figure 2A–D (PDF), the pair alignment table (CSV) and all results (MAT) to `results/`, and prints the key numbers. It takes about 14 min on the test computer, almost all of it for the 100 CCGP draws (`ccgpIterations` at the top of `run_all.m`). `mdsView` sets the viewing angle of the MDS plot.

## Expected output

- 27 of the 200 units are gender-tuned: 26 by the trial-level rule and 1 more by the word-level rule.
- 9 of the 12 pairs are aligned (Benjamini–Hochberg q < 0.05).
- Kruskal stress is 0.178.
- CCGP is 0.602 over 100 draws, with label-shuffle P = 0.218 and unit-shuffle P = 0.248.

These values differ from the published ones because the units are principal components rather than neurons.

## Using the code on other data

Put responses (occurrences × units), `side`, `pair`, `word` and `word_rate` in the same layout as `data/gender_demo.mat`, then run `run_all` or call the functions in `lib/` directly.

## License

MIT; see `LICENSE`.
