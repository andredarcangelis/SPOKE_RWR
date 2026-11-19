#!/usr/bin/env python3
"""
reward.py: score how well a custom SPOKE version reproduces measured radiation responses.

Idea
----
Each radiation schedule (arm) is modelled as a walk on one custom version of SPOKE, started from
seeds that represent that schedule (e.g. the radiation-effect nodes that apply at its dose per
fraction). A good version should give high scores to genes that went UP after that schedule,
low scores to genes that went DOWN, and a stronger response for the schedule with more change.

Ground truth: a CSV with columns arm, gene_symbol, direction (up/down), e.g.
ground_truth/signal2_fig2_genes.csv (SIGNAL 2.0, Fig. 2).

Components (each 0 to 1; 0.5 = no better than chance)
  up_vs_background[arm]  up genes score above other genes (AUROC)
  up_vs_down[arm]        up genes score above down genes (AUROC), when the arm has down genes
  contrast               genes up only in arm B score relatively higher in B's walk than in A's
Reward = mean of the available components.

Hub control: every score is first turned into a percentile WITHIN a degree bin (genes with a
similar number of edges), so a gene doesn't score well just because it has many edges.

Chance check: the ground-truth genes are swapped for random genes of similar degree, many times,
giving a null distribution and an empirical p-value for the reward.

Usage
-----
python reward.py --instance spoke_base \
    --truth ground_truth/signal2_fig2_genes.csv \
    --run 21Gy_1fx=results/spoke_base__arm_1fx.csv \
    --run 30Gy_3fx=results/spoke_base__arm_3fx.csv \
    --contrast 21Gy_1fx,30Gy_3fx [--panel panel_genes.csv] [--n-null 1000]
    [--exclude-groups "epithelial / luminal tumour"]
Run files come from `spokei.py run`. --panel (one column gene_symbol) restricts the background to
the genes actually measured (e.g. the NanoString panel); without it, all SPOKE genes are used.
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import spokei

ROOT = Path(__file__).resolve().parent


def auroc(pos, neg):
    """Probability that a random 'pos' value beats a random 'neg' value (ties count half)."""
    pos, neg = np.asarray(pos, float), np.asarray(neg, float)
    if len(pos) == 0 or len(neg) == 0:
        return np.nan
    allv = np.concatenate([pos, neg])
    ranks = pd.Series(allv).rank().values
    return (ranks[:len(pos)].sum() - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg))


def degree_table(instance):
    cache = ROOT / "results" / f"{instance}__gene_degree.csv"
    if cache.exists():
        return pd.read_csv(cache, dtype={"node_id": str})
    cfg = spokei.load_instance(instance)
    nodes, edges = spokei.build_graph(cfg)
    deg = pd.concat([edges.source, edges.target]).value_counts()
    genes = nodes[nodes.type == "Gene"][["node_id", "name"]].copy()
    genes["degree"] = genes.node_id.map(deg).fillna(0).astype(int)
    cache.parent.mkdir(exist_ok=True)
    genes.to_csv(cache, index=False)
    return genes


def adjusted(scores, bins):
    """Percentile of each gene's score within its degree bin (0 to 1)."""
    return scores.groupby(bins).rank(pct=True)


def components(adj, truth, universe, contrast):
    out = {}
    for arm, g in truth.groupby("arm"):
        up = g[g.direction == "up"].node_id
        down = g[g.direction == "down"].node_id
        a = adj[arm]
        listed = set(truth.node_id)
        bg = universe[~universe.isin(listed)]
        out[f"up_vs_background[{arm}]"] = auroc(a.reindex(up).dropna(), a.reindex(bg).dropna())
        if len(down):
            out[f"up_vs_down[{arm}]"] = auroc(a.reindex(up).dropna(), a.reindex(down).dropna())
    if contrast:
        A, B = contrast
        upA = set(truth[(truth.arm == A) & (truth.direction == "up")].node_id)
        onlyB = [x for x in truth[(truth.arm == B) & (truth.direction == "up")].node_id if x not in upA]
        if onlyB:
            d = (adj[B].reindex(onlyB) - adj[A].reindex(onlyB)).dropna()
            out["contrast"] = float(0.5 + d.mean() / 2) if len(d) else np.nan
    vals = [v for v in out.values() if not np.isnan(v)]
    out["reward"] = float(np.mean(vals)) if vals else np.nan
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--instance", required=True, help="instance used for the runs (for gene degrees)")
    p.add_argument("--truth", required=True)
    p.add_argument("--run", action="append", required=True, help="ARM=path/to/run.csv")
    p.add_argument("--contrast", help="ARM_A,ARM_B: genes up only in B should rise more in B")
    p.add_argument("--panel", help="CSV with gene_symbol column: background = measured genes only")
    p.add_argument("--n-null", type=int, default=1000)
    p.add_argument("--bins", type=int, default=10, help="number of degree bins")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--exclude-groups", default="",
                   help="comma-separated biological_group values to leave out, e.g. "
                        "'epithelial / luminal tumour' (genes that drop because fewer tumour cells were sampled)")
    a = p.parse_args()

    genes = degree_table(a.instance)
    sym2id = genes.drop_duplicates("name").set_index("name").node_id

    truth = pd.read_csv(a.truth)
    truth = truth[truth.direction.isin(["up", "down"])].copy()
    if a.exclude_groups:
        drop = [x.strip() for x in a.exclude_groups.split(",") if x.strip()]
        if "biological_group" not in truth:
            sys.exit("--exclude-groups needs a biological_group column in the ground truth")
        n0 = len(truth); truth = truth[~truth.biological_group.isin(drop)].copy()
        print(f"  excluded {n0 - len(truth)} ground-truth rows in group(s): {'; '.join(drop)}")
    truth["node_id"] = truth.gene_symbol.map(sym2id)
    missing = truth[truth.node_id.isna()]
    if len(missing):
        print(f"  {len(missing)} ground-truth genes not found in this graph, skipped: "
              f"{', '.join(sorted(set(missing.gene_symbol)))}")
    truth = truth.dropna(subset=["node_id"])

    runs = dict(r.split("=", 1) for r in a.run)
    if set(truth.arm) - set(runs):
        sys.exit(f"No --run given for arm(s): {', '.join(sorted(set(truth.arm) - set(runs)))}")

    universe = genes.node_id
    if a.panel:
        panel = pd.read_csv(a.panel).gene_symbol.map(sym2id).dropna()
        universe = universe[universe.isin(set(panel) | set(truth.node_id))]
        print(f"  background restricted to {len(universe):,} panel genes")
    g = genes.set_index("node_id").loc[universe]
    bins = pd.qcut(np.log1p(g.degree), a.bins, labels=False, duplicates="drop")

    adj = {}
    for arm, path in runs.items():
        r = pd.read_csv(path, dtype={"node_id": str}, usecols=["node_id", "score"]).set_index("node_id").score
        adj[arm] = adjusted(r.reindex(universe).fillna(0), bins)

    contrast = tuple(a.contrast.split(",")) if a.contrast else None
    obs = components(adj, truth, pd.Series(universe.values), contrast)

    # chance check: swap each ground-truth gene for a random gene in the same degree bin
    rng = np.random.default_rng(a.seed)
    by_bin = {b: grp.index.values for b, grp in bins.groupby(bins)}
    null = []
    for _ in range(a.n_null):
        t = truth.copy()
        t["node_id"] = [rng.choice(by_bin[bins[x]]) for x in t.node_id]
        null.append(components(adj, t, pd.Series(universe.values), contrast)["reward"])
    null = np.array(null)
    p_val = (1 + np.sum(null >= obs["reward"])) / (1 + len(null))

    print(f"\nInstance: {a.instance}   ground truth: {Path(a.truth).name}   "
          f"genes used: {len(truth)} ({len(truth.gene_symbol.unique())} unique)")
    for k, v in obs.items():
        if k != "reward":
            print(f"  {k:<32} {v:.3f}")
    print(f"  {'REWARD':<32} {obs['reward']:.3f}   (chance: {np.mean(null):.3f}; p = {p_val:.3f}, "
          f"{a.n_null} degree-matched swaps)")

    row = {"instance": a.instance, **obs, "null_mean": float(np.mean(null)), "p_value": float(p_val)}
    out = ROOT / "results" / "reward_table.csv"
    pd.DataFrame([row]).to_csv(out, mode="a", header=not out.exists(), index=False)
    print(f"  appended to {out.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
