# spoke-instances: custom SPOKE instances

Build variants of SPOKE, run the same seeds on each, and compare the rankings.
This folder sits inside the NASA SPOKE folder; the rules in ../CLAUDE.md
(read-only Pitt, tagging, citations) still apply.

## Concepts
- **Base**: plain SPOKE as CSV in `../spoke_export/` (`nodes.csv`: node_id,name,type;
  `edges.csv`: source,target,edge_type). Never edit these files.
  Node IDs there are `<Type>:<SPOKE identifier>`, e.g. `Anatomy:UBERON:0003233`.
- **Overlay**: a set of added nodes/edges in `overlays/<name>/nodes.csv` and `edges.csv`.
  Edge columns: source,target,edge_type,weight,citation. New node IDs start with `ext:`.
  Edges to SPOKE nodes must use the exact node_id from `../spoke_export/nodes.csv`
  (look it up there; don't build it by hand from a name).
- **Instance**: `instances/<name>.toml` = base + list of overlays + optional
  `edge_type_weights` (multipliers) and `drop_edge_types`. Copy an existing file to start one.
- **Seeds**: `seeds/<name>.csv` with node_id,weight (e.g. log2 fold changes).

## Commands (run from this folder, in the spoke_env environment)
- `python spokei.py list` lists the instances
- `python spokei.py check <instance>` validates it; catches edges pointing at IDs that don't exist
- `python spokei.py run <instance> --seeds seeds/<file>.csv` personalized PageRank, saved to results/
- `python spokei.py compare results/<a>.csv results/<b>.csv` Spearman, cosine, top-k overlap per type, biggest movers

## How Claude should handle requests
- "Make an instance that ..." : write the .toml; create the overlay folder if needed.
- "Add X to overlay Y" : build it in spoke-write first (tag `added_by: "fractionation_ext"` and
  `overlay: "<name>"`), then export: query spoke-write for nodes/edges with that overlay tag and
  write them to overlays/<name>/ in the format above. Copied SPOKE nodes are not re-exported
  (they're already in the base); only the new nodes and all new edges are.
- Store the citation in the Neo4j property `citation` (the CSV column `source` is taken by the edge endpoint).
- Always run `check` before `run`, and `run` the base instance with the same seeds before comparing.
- Report results in plain language: what moved, by how much, and whether it's more than noise.
  A single run is not a significance test; say so when it matters.
- The demo_* instances and data/toy/ are a small fake graph for testing the tool only.
  Never present demo results as SPOKE results.

## Reward: does a version reproduce measured radiation responses?
- `ground_truth/signal2_fig2_genes.csv`: genes labelled in SIGNAL 2.0 Fig. 2 (21 Gy x 1 vs 10 Gy x 3,
  breast, paired biopsies). Only the labelled subset of the DEGs; log2fc_approx is read off the plot.
  LBP and MARCO are text-only (arm/direction unknown) and are skipped by the reward.
- `python reward.py --instance <inst> --truth ground_truth/signal2_fig2_genes.csv
   --run 21Gy_1fx=results/<inst>__<arm1 seeds>.csv --run 30Gy_3fx=results/<inst>__<arm3 seeds>.csv
   --contrast 21Gy_1fx,30Gy_3fx`
  Each arm's run uses seeds that represent that schedule (the radiation-effect nodes that apply
  at its dose). Scores are degree-adjusted; the p-value comes from degree-matched gene swaps.
- Compare versions only on the SAME seeds and ground truth. Never use SIGNAL 2.0 as a source of
  edges or effect nodes; it is the test set.
- `ground_truth/toy_planted.csv` / `toy_scrambled.csv` + `seeds/toy_arm*.csv` on demo_base are the
  self-test: planted should give a high reward with a small p, scrambled about 0.5.
