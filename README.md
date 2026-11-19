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

The base SPOKE export (`nodes.csv` and `edges.csv`) is not in this repository. It is ~1 GB, above GitHub's file-size limit, and SPOKE data should not be redistributed without permission from the SPOKE team. Anyone with SPOKE access can make their own export (step 4 below).

`data/toy/` and the `demo_*` instances are a small fake graph for testing the tool only.

## Setup with Claude Code and MCP

### 1. Folder layout

`spokei.py` expects the export in a sibling folder:

```
NASA SPOKE/
├── CLAUDE.md          # copy of CLAUDE_for_NASA_SPOKE_folder.md (project rules)
├── spoke_export/      # nodes.csv, edges.csv (step 4)
└── SPOKE_RWR/         # this repository
```

```
git clone https://github.com/andredarcangelis/SPOKE_RWR.git
cp SPOKE_RWR/CLAUDE_for_NASA_SPOKE_folder.md CLAUDE.md
```

### 2. Python

Python 3.11+ (uses `tomllib`), plus numpy, pandas and scipy:

```
python -m venv spoke_env && source spoke_env/bin/activate
pip install numpy pandas scipy
```

### 3. Connect Claude Code to SPOKE

You need access to a SPOKE Neo4j server (ask the SPOKE team or whoever manages your group's server). If it is reached through an SSH tunnel, open the tunnel first, mapping the server's Bolt port to local port 17687:

```
ssh -N -L 17687:localhost:7687 <user>@<spoke-server>
```

Then add a read-only Neo4j MCP connection to Claude Code (needs `uv` installed). Run this from the `NASA SPOKE` folder:

```
claude mcp add spoke-pitt \
  -e NEO4J_URI=bolt://localhost:17687 \
  -e NEO4J_USERNAME=<username> \
  -e NEO4J_PASSWORD=<password> \
  -e NEO4J_DATABASE=spoke \
  -e NEO4J_READ_ONLY=true \
  -- uvx mcp-neo4j-cypher@0.6.0 --transport stdio
```

Restart Claude Code, then check with `/mcp`. Claude can now look up SPOKE nodes and identifiers. The connection name `spoke-pitt` and database `spoke` match the project's CLAUDE.md; change them if your server differs.

Optional: to build overlays, add a second connection named `spoke-write` to a local Neo4j database (for example in Neo4j Desktop on `bolt://localhost:7687`), without `NEO4J_READ_ONLY`. CLAUDE.md explains how overlays are built there and exported to `overlays/<name>/`.

### 4. Export SPOKE to CSV

The MCP connection is for lookups; a full export (millions of edges) should be done with a script using the Neo4j Python driver, not through MCP. You can ask Claude Code to write and run it with a prompt like:

```
Read-only on spoke-pitt. Write and run a Python script (neo4j driver) that exports SPOKE to
../spoke_export/nodes.csv (node_id,name,type) and ../spoke_export/edges.csv (source,target,edge_type).
- node_id = "<Label>:<identifier>", e.g. Anatomy:UBERON:0003233
- Keep only human proteins (Protein nodes from Homo sapiens).
- Node labels: <list the labels you want, or "all">.
- Keep duplicate edge rows; record the labels, filters and counts in ../spoke_export/export_settings.txt.
- Stream in batches; report node and edge counts per label and edge type when done.
```

For reference, the human-protein export of all 28 node types used in this project (October 2026) had 3,684,631 nodes and 12,592,187 edges. The 17-type version had 3,145,112 nodes.

### 5. Check the tool works

```
cd SPOKE_RWR
python spokei.py check demo_base
python spokei.py run demo_base --seeds seeds/toy_arm1.csv
python spokei.py check spoke_base     # needs the export from step 4
```

`ground_truth/toy_planted.csv` and `toy_scrambled.csv` with `seeds/toy_arm*.csv` on `demo_base` are the self-test for `reward.py`: planted should score high with a small p-value, scrambled about 0.5.

Run results on full SPOKE are ~480 MB each, so avoid committing them.
