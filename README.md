# spoke-instances

Build custom versions ("instances") of the SPOKE knowledge graph, run the same random walk with restart on each, and compare the rankings.

An instance is the base SPOKE export plus optional overlays (added nodes and edges), edge-type weight multipliers, or removed edge types. Each instance is one small TOML file in `instances/`.

## The random walk (`spokei.py`)

- Personalised PageRank (random walk with restart), solved by power iteration until the total change is below 1e-10 (at most 1,000 iterations).
- Restart probability 0.1 by default (`restart` in the instance file).
- Edges are undirected by default (`directed = false`).
- Duplicate edge rows are summed as weight unless `collapse_duplicate_edges = true`.
- Mass at dangling nodes is returned to the seeds.
- Optional `hub_penalty` (default 0): scales each edge into a node by degree^(-hub_penalty), making high-degree nodes less likely targets.
- Seeds can carry weights (e.g. log2 fold change), including negative ones; because the walk is linear, the result is the weighted sum of per-seed walks.

## Commands

```
python spokei.py list                                   # show instances
python spokei.py check   <instance>                     # validate, print graph size
python spokei.py run     <instance> --seeds <file.csv>  # random walk -> results/
python spokei.py compare <result_a.csv> <result_b.csv>  # how two rankings differ
```

`reward.py` scores a ranking against measured radiation-response genes (see `ground_truth/`).

## File formats (CSV, header required)

- Base and overlay `nodes.csv`: `node_id, name, type`
- Base and overlay `edges.csv`: `source, target, edge_type[, weight]` (weight defaults to 1)
- `seeds/*.csv`: `node_id, weight`

## Data not included

The base SPOKE export (`../spoke_export/nodes.csv` and `edges.csv`) is not in this repository. It comes from a SPOKE Neo4j instance; you need your own export in the format above.

`data/toy/` and the `demo_*` instances are a small fake graph for testing the tool only.

## Requirements

Python 3.11+ (uses `tomllib`), numpy, pandas, scipy.
