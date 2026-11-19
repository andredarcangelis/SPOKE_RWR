# NASA SPOKE: radiation fractionation extension

## What this project is
Adding radiation fractionation to the SPOKE knowledge graph. Today, a tumor treated with
2 Gy x 30 and one treated with 8 Gy x 3 would both just link to one "radiotherapy" node.
The goal is fractionation-regime nodes connected to existing SPOKE nodes (genes, GO terms),
so a schedule changes *which* biology a tumor connects to, not only how strongly.

## Databases (two MCP connections)
- **spoke-pitt**: Pitt's shared SPOKE graph (database `spoke`), via SSH tunnel on
  localhost:17687. READ ONLY. Use it to look up SPOKE nodes and their identifiers.
  Never try to write here. APOC is not installed, so the schema tool will fail; use plain queries.
- **spoke-write**: local Neo4j Desktop server "Radiation Fractionation" (database `neo4j`)
  on localhost:7687. This is the only place to create nodes and edges.

If spoke-pitt is not connected, the tunnel is probably closed; tell me rather than guessing.

## Rules for every write to spoke-write
1. Tag every node and edge you create with `added_by: "fractionation_ext"`.
2. Every new edge needs a `citation` property: a PMID, DOI, or dataset ID (e.g. OSD-xxx, GSExxxxx).
   If there is no citation, do not create the edge; list it for me instead.
3. When copying a node from Pitt, keep its original label and identifier properties
   unchanged, and add `copied_from: "pitt_spoke"`, so it can be matched back to SPOKE.
4. Use MERGE on the identifier, not CREATE, so running a step twice does not make duplicates.
5. Never delete untagged data. To undo our additions only:
   `MATCH (n {added_by: "fractionation_ext"}) DETACH DELETE n`
6. Before any write, say in one or two plain sentences what it will add, then run it.

## Regime definitions
Source: ASTRO/ASCO/AUA guideline, Hypofractionated Radiation Therapy for Localized Prostate
Cancer (2018). Use this as the `citation` on the is-a edges.
- Conventional: 1.8 to 2.0 Gy per fraction (180 to 200 cGy)
- Moderate hypofractionation: 2.4 to 3.4 Gy per fraction (240 to 340 cGy)
- Ultrahypofractionation: >= 5 Gy per fraction (>= 500 cGy)
These are prostate-specific, and the guideline leaves gaps (2.0 to 2.4 Gy and 3.4 to 5 Gy are
not assigned). Do not silently bin a schedule that falls in a gap; flag it for me.
Each regime node links to a general Radiotherapy node with an is-a edge.

## Edge options, in order of strength
1. Label only: regime node linked to the tumors on that schedule (no biology edges).
2. Literature: regime node linked to SPOKE genes/GO terms, one citation per edge (curated by me).
3. Data: regime node linked to genes that change in studies using that regime
   (differential expression from GeneLab or GEO).

## How to talk to me
I'm a biochemistry student; explain things plainly and directly. When you run a query,
summarize what it found in plain language. Say plainly when something could not be confirmed.
Do not invent citations, gene–regime links, or dose thresholds.
