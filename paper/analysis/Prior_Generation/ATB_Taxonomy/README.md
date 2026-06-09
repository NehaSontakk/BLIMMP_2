Resolves species names for all ATB genomes that have been tantan-masked into full
NCBI taxonomic lineages. Produces the per-sample lineage table consumed by
`lineage_specific_ko_matrix.py` downstream.

---

## Inputs

| Source | Path (in script) | Description |
|--------|-----------------|-------------|
| Tantan-masked genomes | `ATB_TANTAN/*_masked.fa` | Glob of masked FASTA files; used only to extract sample IDs — file contents are not read |
| Master file list | `DOWNLOAD_SCRIPTS/file_list.all.20240805.tsv` | Two columns used: `sample` and `species_sylph` (sylph-assigned species name, GTDB-derived) |
| NCBI taxonomy dump | `ATB_Taxonomy/` (taxonkit `--data-dir`) | Standard taxonkit data directory (`names.dmp`, `nodes.dmp`, etc.) — must be pre-downloaded |

---

## Outputs

All files are written to the working directory (wherever the script is run from).
All TSVs are written **without a header row**.

| File | Columns | Description |
|------|---------|-------------|
| `ATB_species_clean_map.tsv` | `species_sylph`, `species_clean` | Audit map of raw sylph name → cleaned name; deduplicated |
| `ATB_sample_species.tsv` | `sample`, `species_clean` | One row per genome with cleaned species name |
| `ATB_sample_species1.tsv` | `sample`, `species_clean`, `taxid` | As above plus resolved NCBI taxid; this is the direct input to `taxonkit lineage` |
| `ATB_sample_species_full_lineage.tsv` | *(raw taxonkit output)* | Intermediate — full lineage string + ranks, may contain malformed "cellular root" lines |
| `ATB_sample_species_full_lineage_clean.tsv` | same, cleaned | **Final output** — "cellular root" header lines stripped; consumed downstream |

---

## What it does (step by step)

### 1. Identify the working genome set
Globs `ATB_TANTAN/*_masked.fa` and strips `_masked.fa` to get sample IDs. Only
samples present in this directory are carried forward — this keeps the taxonomy
file in sync with genomes that have actually been processed through tantan masking.

### 2. Filter the master list
Loads only the `sample` and `species_sylph` columns from the master TSV and
keeps rows whose `sample` appears in the tantan set.

### 3. Clean species names (`clean_taxon_name`)
`sylph` reports GTDB-style names which do not map directly to NCBI taxonomy.
Three transformations are applied:

| Pattern | Example in → out | Reason |
|---------|-----------------|--------|
| Genus partition suffix (`_A`, `_B`, `_K`, `_G3`, …) | `Bacillus_A` → `Bacillus` | GTDB splits clades with letter suffixes; NCBI has none |
| Species epithet partition suffix | `thuringiensis_S` → `thuringiensis` | Same GTDB convention |
| Numeric sp. placeholder (`sp` + 3+ digits) | `sp010998615` → `sp.` | sylph emits GTDB accession-based placeholders for unnamed species |
| `Candidatus` / `Ca.` prefix | preserved; genus normalization starts at next token | These are valid NCBI names and must not be stripped |

The raw→clean mapping is written to `ATB_species_clean_map.tsv` for auditing.

### 4. Bulk taxid resolution (`name2taxids_bulk`)
All **unique** cleaned names are sent to `taxonkit name2taxid` in a single
subprocess call (one name per line). This is the primary performance optimization
over the original script, which called taxonkit once per genome row.

Output is a dict `name → [taxid, taxid2, …]`. A name can return multiple taxids
when it is ambiguous in NCBI (e.g. a genus name that maps to both a bacterial and
an archaeal node).

### 5. LCA collapse (`lca_bulk`)
For names with multiple taxid candidates, `taxonkit lca` is called in bulk to
compute the Lowest Common Ancestor — the most specific unambiguous NCBI node.
Singletons pass through unchanged. Empty lists (lookup failures) produce an
empty taxid string.

> **No fallback for failed lookups.** If a name fails `name2taxid` entirely
> (e.g. a genus+sp. combination that doesn't exist in NCBI under any form),
> the taxid is recorded as `""`. These rows will produce an empty lineage in the
> final output. Check the rows with empty taxids in
> `ATB_sample_species1.tsv` if lineage coverage looks low.

### 6. Write intermediate TSVs
`ATB_sample_species.tsv` (sample + clean name) and `ATB_sample_species1.tsv`
(sample + clean name + taxid) are written. The latter is the direct input to
`taxonkit lineage`.

### 7. Run `taxonkit lineage`
Called as a shell subprocess:
```
taxonkit lineage \
  --data-dir <TAXONKIT_DATA> \
  --taxid-field 3 \
  --show-lineage-ranks \
  -j <THREADS> \
  ATB_sample_species1.tsv \
  > ATB_sample_species_full_lineage.tsv
```
`--taxid-field 3` tells taxonkit that column 3 (1-indexed) contains the taxid.
`--show-lineage-ranks` adds a parallel column naming the rank of each node in
the lineage string (e.g. `superkingdom;phylum;class;…`).

### 8. Strip malformed lines
Lines starting with `"cellular root"` are artifacts of taxonkit when a taxid
resolves to the root node (taxid 1). These are removed to produce
`ATB_sample_species_full_lineage_clean.tsv`.


## Downstream
`ATB_sample_species_full_lineage_clean.tsv` is consumed by
`lineage_specific_ko_matrix.py` in `ATB_KO_Frequencies/` to group genomes by
taxonomic lineage before computing per-lineage KO frequency priors.
