## Reward: does a version reproduce measured radiation responses?
- `ground_truth/signal2_fig2_genes.csv`: genes labelled in SIGNAL 2.0 Fig. 2 (21 Gy x 1 vs 10 Gy x 3,
  breast, paired biopsies). Only the labelled subset of the DEGs; log2fc_approx is read off the plot.
  LBP and MARCO are text-only (arm/direction unknown) and are skipped by the reward.
- `python reward.py --instance <inst> --truth ground_truth/signal2_fig2_genes.csv
   --run 21Gy_1fx=results/<inst>__<arm1 seeds>.csv --run 30Gy_3fx=results/<inst>__<arm3 seeds>.csv
   --contrast 21Gy_1fx,30Gy_3fx`
  Add `--exclude-groups "epithelial / luminal tumour"` to leave out genes that likely dropped only
  because the post-RT lumpectomy samples contain fewer tumour cells; report both versions.
  Each arm's run uses seeds that represent that schedule (the radiation-effect nodes that apply
  at its dose). Scores are degree-adjusted; the p-value comes from degree-matched gene swaps.
- Compare versions only on the SAME seeds and ground truth. Never use SIGNAL 2.0 as a source of
  edges or effect nodes; it is the test set.
- `ground_truth/toy_planted.csv` / `toy_scrambled.csv` + `seeds/toy_arm*.csv` on demo_base are the
  self-test: planted should give a high reward with a small p, scrambled about 0.5.
