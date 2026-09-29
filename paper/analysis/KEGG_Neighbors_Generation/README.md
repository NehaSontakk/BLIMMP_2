# KEGG Neighbors Generation

Post-processing scripts that reformat HPC-generated neighbor and reaction data for downstream consumers.

The core one-hop/two-hop/module-all-hop neighbor statistics (`MODULE_ALL_NEIGHBOR_DATA/`, `ONE_HOP_NEIGHBOR_DATA-2/`, `TWO_HOP_NEIGHBOR_DATA-2/`) are generated on the HPC — see `Prior_Generation/`, not this folder.

**Contents:**
- `KO_to_Reaction_Step_Mapping.ipynb` — builds `KO_to_Reactions.xlsx`, the KO-to-reaction mapping, from `KEGG_Prokaryote_Modules_with_Definitions.xlsx` (see `KEGG_Graph_Generation/`)
- `Generating_reverse_conditional_json_for_BLIMMP-Web.ipynb` — reformats `MODULE_ALL_NEIGHBOR_DATA/Module_AllHop_Refilled_domain_level_priors.json` into the conditional-probability JSON format BLIMMP-Web consumes
- `KO_to_Reactions.xlsx` — output of the mapping notebook above


## Neighbor data archive

`KEGG_Neighbor_Data.zip` contains the actual one-hop, two-hop, and module-all-hop neighbor/co-occurrence data BLIMMP consumes (generated on the HPC via `Prior_Generation/`, not by the notebooks in this folder): `MODULE_ALL_NEIGHBOR_DATA/`, `ONE_HOP_NEIGHBOR_DATA-2/`, `TWO_HOP_NEIGHBOR_DATA-2/`, `Module_Neighborhood/`, `One_Hop_Neighbor.txt`, `Two_Hop_Neighbor.txt`. The live copies BLIMMP actually reads from during pipeline runs are not this archive — see `blimmp_src/BLIMMP_8Sep2026/BLIMMP/`.
