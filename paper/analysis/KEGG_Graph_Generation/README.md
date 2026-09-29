# KEGG Graph Generation

Parses KEGG metabolic module definitions into the directed acyclic graphs BLIMMP uses for step and module inference.

**Pipeline order:**
1. `Fetch_Prokaryote_Modules_From_KEGG.ipynb` — retrieves the list of prokaryote-relevant KEGG modules
2. `Fetch_Definitions_From_KEGG.ipynb` — retrieves each module's logical definition string
3. `GRAPH_GENERATION_USING_KEGG_CFG_Updated_Aug25.ipynb` — parses each definition with the grammar in `new_updated_grammar.txt` and builds per-module node/adjacency graphs
4. `Aggregate_Module_Info_Updated_13AUG2026.ipynb` — merges the generated graphs with AllTheBacteria domain-level KO frequencies and edge priors into the node/adjacency JSONs used by BLIMMP-Explorer

**Inputs:** KEGG REST API (module list and definitions fetched once, not re-pulled at runtime), `new_updated_grammar.txt`, AllTheBacteria domain-level priors (see `Prior_Generation/`).

**Outputs:** per-module `module_*_nodes.json` / `module_*_adjacency.json` graphs (consumed directly by the BLIMMP pipeline) and `all_module_nodes_13AUG2026.json` / `all_module_adjacency_links_13AUG2026.json` (consumed by BLIMMP-Explorer).


## Generated graphs archive

`KEGG_Graphs_Generated_13AUG2026.tar.gz` is a compressed snapshot of all 1361 generated files (per-module `nodes.json`, `adjacency.json`, `paths.json`, and `graph.png`), for reference and reproducibility checking. Decompress with `tar -xzf KEGG_Graphs_Generated_13AUG2026.tar.gz` to reproduce the exact directory structure BLIMMP's `module_json_dir` config expects. The live copy BLIMMP actually reads from during pipeline runs is not this archive — see `blimmp_src/BLIMMP_8Sep2026/BLIMMP/KEGG_Graphs_Generated_13AUG2026/`.


## Generated graphs archive

`KEGG_Graphs_Generated_13AUG2026.tar.gz` is a compressed snapshot of all 1361 generated files (per-module `nodes.json`, `adjacency.json`, `paths.json`, and `graph.png`), for reference and reproducibility checking. Decompress with `tar -xzf KEGG_Graphs_Generated_13AUG2026.tar.gz` to reproduce the exact directory structure BLIMMP's `module_json_dir` config expects. The live copy BLIMMP actually reads from during pipeline runs is not this archive — see `blimmp_src/BLIMMP_8Sep2026/BLIMMP/KEGG_Graphs_Generated_13AUG2026/`.


## Generated graphs archive

`KEGG_Graphs_Generated_13AUG2026.tar.gz` is a compressed snapshot of all 1361 generated files (per-module `nodes.json`, `adjacency.json`, `paths.json`, and `graph.png`), for reference and reproducibility checking. Decompress with `tar -xzf KEGG_Graphs_Generated_13AUG2026.tar.gz` to reproduce the exact directory structure BLIMMP's `module_json_dir` config expects. The live copy BLIMMP actually reads from during pipeline runs is not this archive — see `blimmp_src/BLIMMP_8Sep2026/BLIMMP/KEGG_Graphs_Generated_13AUG2026/`.


## Generated graphs archive

`KEGG_Graphs_Generated_13AUG2026.tar.gz` is a compressed snapshot of all 1361 generated files (per-module `nodes.json`, `adjacency.json`, `paths.json`, and `graph.png`), for reference and reproducibility checking. Decompress with `tar -xzf KEGG_Graphs_Generated_13AUG2026.tar.gz` to reproduce the exact directory structure BLIMMP's `module_json_dir` config expects. The live copy BLIMMP actually reads from during pipeline runs is not this archive — see `blimmp_src/BLIMMP_8Sep2026/BLIMMP/KEGG_Graphs_Generated_13AUG2026/`.
