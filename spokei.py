#!/usr/bin/env python3
"""
spokei: build, run and compare custom SPOKE instances.

An instance = the base SPOKE graph + zero or more overlays (your added nodes/edges)
+ optional edge-type weights or removals. Each instance is one small TOML file in
instances/. Nothing here touches Neo4j; Claude (via the spoke-write MCP) builds
overlays in Neo4j and exports them to overlays/<name>/ as CSV.

Commands
  python spokei.py list                                   show instances
  python spokei.py check   <instance>                     validate it, print graph size
  python spokei.py run     <instance> --seeds <file.csv>  personalized PageRank -> results/
  python spokei.py compare <result_a.csv> <result_b.csv>  how rankings differ

File formats (CSV, header required)
  base / overlay nodes.csv : node_id, name, type
  base / overlay edges.csv : source, target, edge_type[, weight]   (weight defaults to 1)
  seeds/*.csv              : node_id, weight   (e.g. log2 fold change; may be negative)
"""

import argparse
import sys
import tomllib
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

ROOT = Path(__file__).resolve().parent
INSTANCES, OVERLAYS, RESULTS = ROOT / "instances", ROOT / "overlays", ROOT / "results"
DEFAULTS = {"base": "../spoke_export", "overlays": [], "directed": False, "restart": 0.1,
            "edge_type_weights": {}, "drop_edge_types": [], "description": "",
            "collapse_duplicate_edges": False,
            "hub_penalty": 0.0}


# ---------------------------------------------------------------- loading
def load_instance(name):
    path = INSTANCES / f"{name}.toml"
    if not path.exists():
        sys.exit(f"No instance called '{name}' (looked for {path}).")
    with open(path, "rb") as f:
        cfg = {**DEFAULTS, **tomllib.load(f)}
    cfg["name"] = name
    return cfg


def read_csv(path, required):
    if not Path(path).exists():
        sys.exit(f"File not found: {path}")
    df = pd.read_csv(path, dtype={c: str for c in ("node_id", "source", "target")})
    missing = [c for c in required if c not in df.columns]
    if missing:
        sys.exit(f"{path} is missing column(s): {', '.join(missing)}")
    return df


def build_graph(cfg):
    """Return (nodes DataFrame, edges DataFrame) for an instance, with checks."""
    base = ROOT / cfg["base"]
    nodes = [read_csv(base / "nodes.csv", ["node_id", "name", "type"]).assign(origin="base")]
    edges = [read_csv(base / "edges.csv", ["source", "target", "edge_type"]).assign(origin="base")]
    for ov in cfg["overlays"]:
        d = OVERLAYS / ov
        if not d.exists():
            sys.exit(f"Overlay '{ov}' not found at {d}")
        if (d / "nodes.csv").exists():
            nodes.append(read_csv(d / "nodes.csv", ["node_id", "name", "type"]).assign(origin=ov))
        if (d / "edges.csv").exists():
            edges.append(read_csv(d / "edges.csv", ["source", "target", "edge_type"]).assign(origin=ov))
    nodes = pd.concat(nodes, ignore_index=True)
    edges = pd.concat(edges, ignore_index=True)
    if cfg["collapse_duplicate_edges"]:
        # identical (source, target, edge_type) rows count once instead of being summed
        before = len(edges)
        edges = edges.drop_duplicates(subset=["source", "target", "edge_type"]).reset_index(drop=True)
        print(f"  collapsed {before - len(edges):,} duplicate edge rows")

    dup = nodes[nodes.duplicated("node_id", keep=False)]
    if len(dup):
        sys.exit(f"Overlay node IDs clash with existing nodes:\n{dup.head(10).to_string(index=False)}")

    known = set(nodes.node_id)
    bad = edges[~edges.source.isin(known) | ~edges.target.isin(known)]
    if len(bad):
        sys.exit("These edges point to node IDs that are not in the base or overlay nodes "
                 f"({len(bad)} total; first 10 shown). Check the IDs match the base nodes.csv:\n"
                 f"{bad.head(10)[['source', 'target', 'edge_type', 'origin']].to_string(index=False)}")

    if "weight" not in edges:
        edges["weight"] = 1.0
    edges["weight"] = edges["weight"].fillna(1.0).astype(float)
    if (edges.weight < 0).any():
        sys.exit("Edge weights must be >= 0 (a random walk can't follow a negative edge).")
    edges = edges[~edges.edge_type.isin(cfg["drop_edge_types"])].copy()
    for et, w in cfg["edge_type_weights"].items():
        edges.loc[edges.edge_type == et, "weight"] *= float(w)
    return nodes, edges


# ---------------------------------------------------------------- PageRank
def personalized_pagerank(nodes, edges, seeds, restart, directed, tol=1e-10, max_iter=1000, hub_penalty=0.0):
    """Random walk with restart. Each step: with prob `restart` jump to the seeds
    (in proportion to |seed weight|, signed), otherwise follow an edge in proportion to its weight.
    Because PPR is linear, signed seeds give the FC-weighted sum of per-gene walks."""
    idx = pd.Series(np.arange(len(nodes)), index=nodes.node_id.values)
    s, t, w = idx[edges.source].values, idx[edges.target].values, edges.weight.values
    if not directed:
        s, t, w = np.concatenate([s, t]), np.concatenate([t, s]), np.concatenate([w, w])
    n = len(nodes)
    A = sp.csr_matrix((w, (s, t)), shape=(n, n))           # duplicates are summed
    if hub_penalty:
        # make stepping INTO a node less likely the more edges it has: weight x degree^(-hub_penalty)
        deg = np.bincount(np.concatenate([s, t]) if directed else s, minlength=n).astype(float)
        scale = np.where(deg > 0, deg, 1.0) ** (-hub_penalty)
        A = (A @ sp.diags(scale)).tocsr()
    out = np.asarray(A.sum(axis=1)).ravel()
    dangling = out == 0
    inv = np.where(dangling, 0.0, 1.0 / np.where(dangling, 1, out))
    PT = (sp.diags(inv) @ A).T.tocsr()                     # transpose of row-stochastic P

    seed_vec = np.zeros(n)
    missing = seeds[~seeds.node_id.isin(idx.index)]
    if len(missing):
        print(f"  warning: {len(missing)} seed node(s) not in this graph, skipped "
              f"(e.g. {', '.join(missing.node_id.head(5))})")
    seeds = seeds[seeds.node_id.isin(idx.index)]
    np.add.at(seed_vec, idx[seeds.node_id].values, seeds.weight.values.astype(float))
    if np.abs(seed_vec).sum() == 0:
        sys.exit("No usable seeds.")
    seed_vec /= np.abs(seed_vec).sum()

    p = seed_vec.copy()
    for it in range(max_iter):
        new = (1 - restart) * (PT @ p + seed_vec * p[dangling].sum()) + restart * seed_vec
        if np.abs(new - p).sum() < tol:
            return new, it + 1
        p = new
    print(f"  warning: did not converge in {max_iter} iterations")
    return p, max_iter


# ---------------------------------------------------------------- commands
def cmd_list(_):
    files = sorted(INSTANCES.glob("*.toml"))
    if not files:
        print("No instances yet. Add a .toml file to instances/.")
    for f in files:
        cfg = load_instance(f.stem)
        ov = ", ".join(cfg["overlays"]) or "none"
        print(f"{f.stem:<22} overlays: {ov:<30} {cfg['description']}")


def cmd_check(a):
    cfg = load_instance(a.instance)
    nodes, edges = build_graph(cfg)
    print(f"Instance '{a.instance}' is valid.")
    print(f"  nodes: {len(nodes):,}   edges: {len(edges):,}   directed: {cfg['directed']}   restart: {cfg['restart']}")
    added = edges[edges.origin != "base"]
    if len(added):
        print("  added edges by overlay and type:")
        print(added.groupby(["origin", "edge_type"]).size().rename("n").to_string())
    print(f"  duplicate edge rows collapsed: {cfg['collapse_duplicate_edges']}   hub penalty: {cfg['hub_penalty']}")
    if cfg["edge_type_weights"] or cfg["drop_edge_types"]:
        print(f"  reweighted types: {cfg['edge_type_weights'] or 'none'}   dropped types: {cfg['drop_edge_types'] or 'none'}")


def cmd_run(a):
    cfg = load_instance(a.instance)
    nodes, edges = build_graph(cfg)
    seeds = read_csv(a.seeds, ["node_id", "weight"])
    print(f"Running '{a.instance}' with {len(seeds)} seeds from {a.seeds} ...")
    score, iters = personalized_pagerank(nodes, edges, seeds, cfg["restart"], cfg["directed"], hub_penalty=cfg["hub_penalty"])
    res = nodes[["node_id", "name", "type", "origin"]].assign(score=score)
    res["rank_in_type"] = res.groupby("type")["score"].rank(ascending=False, method="min").astype(int)
    res = res.sort_values("score", ascending=False)
    RESULTS.mkdir(exist_ok=True)
    out = RESULTS / f"{a.instance}__{Path(a.seeds).stem}.csv"
    res.to_csv(out, index=False)
    print(f"  converged in {iters} iterations -> {out.relative_to(ROOT)}")
    print(res.head(a.top)[["name", "type", "score", "rank_in_type"]].to_string(index=False))


def cmd_compare(a):
    A, B = pd.read_csv(a.a, dtype={"node_id": str}), pd.read_csv(a.b, dtype={"node_id": str})
    m = A.merge(B, on="node_id", suffixes=("_a", "_b"))          # nodes present in both
    only_b = B[~B.node_id.isin(A.node_id)]
    print(f"A = {Path(a.a).name}\nB = {Path(a.b).name}")
    print(f"Shared nodes: {len(m):,}   only in B (added by its overlays): {len(only_b):,}\n")

    va, vb = m.score_a.values, m.score_b.values
    cos = va @ vb / (np.linalg.norm(va) * np.linalg.norm(vb))
    rho = m.score_a.rank().corr(m.score_b.rank())               # Spearman
    print(f"Whole graph:  Spearman rank correlation {rho:.4f}   cosine similarity {cos:.6f}")
    print("  (1.0 = identical. Cosine near 1 with lower Spearman = same overall shape, reshuffled order.)\n")

    print(f"Top-{a.top} overlap by node type (shared nodes only):")
    for t, g in m.groupby("type_a"):
        ta = set(g.nlargest(a.top, "score_a").node_id)
        tb = set(g.nlargest(a.top, "score_b").node_id)
        print(f"  {t:<22} {len(ta & tb):>3}/{a.top}")

    m["rank_change"] = m.rank_in_type_a - m.rank_in_type_b         # positive = moved up in B
    movers = m.reindex(m.rank_change.abs().sort_values(ascending=False).index).head(a.top)
    print(f"\nBiggest rank movers within their type (positive = higher in B):")
    print(movers[["name_a", "type_a", "rank_in_type_a", "rank_in_type_b", "rank_change"]]
          .rename(columns={"name_a": "name", "type_a": "type", "rank_in_type_a": "rank_A",
                           "rank_in_type_b": "rank_B"}).to_string(index=False))
    out = RESULTS / f"compare__{Path(a.a).stem}__vs__{Path(a.b).stem}.csv"
    m.to_csv(out, index=False)
    print(f"\nFull comparison table -> {out.relative_to(ROOT)}")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    sub.add_parser("list").set_defaults(f=cmd_list)
    c = sub.add_parser("check"); c.add_argument("instance"); c.set_defaults(f=cmd_check)
    r = sub.add_parser("run"); r.add_argument("instance"); r.add_argument("--seeds", required=True)
    r.add_argument("--top", type=int, default=15); r.set_defaults(f=cmd_run)
    k = sub.add_parser("compare"); k.add_argument("a"); k.add_argument("b")
    k.add_argument("--top", type=int, default=10); k.set_defaults(f=cmd_compare)
    a = p.parse_args()
    a.f(a)


if __name__ == "__main__":
    main()
