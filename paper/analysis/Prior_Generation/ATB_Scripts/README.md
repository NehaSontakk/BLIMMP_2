Builds a binary KO × genome presence/absence matrix across all ~21,000 ATB
genomes, plus two global summary tables (per-KO counts and KO co-occurrence).
These are the foundation for all downstream lineage-specific frequency priors.

## Scripts

| Script | Input | Output |
|--------|-------|--------|
| `create_ko_matrix.py` | `PROCESS_HMMER_OUTPUT/*_HMMER_grouped_hits_dedup.csv` | `ko_matrix.tsv` |
| `create_ko_counts.py` | `ko_matrix.tsv` | `ko_counts.tsv` |
| `create_ko_intersection.py` | `ko_matrix.tsv` | `ko_intersection.tsv` |
| `run_ko_matrix.sh` | — | SLURM wrapper for `create_ko_matrix.py` |

---

## `create_ko_matrix.py`

Pass 1 reads all input CSVs once to collect the union of every
KO ID seen across all genomes (no single genome defines the full KO space).
Pass 2 reads each file again and marks a 1 for any KO with ≥1 hit below the
E-value threshold (`1e-5`). Copy number is ignored, presence/absence only.

**Output format:** rows = KO IDs, columns = sample IDs, values ∈ {0, 1}, `uint8`.
Header row and named index column (`KO id`) are included.


- Sample ID is the filename with `_HMMER_grouped_hits_dedup.csv` stripped. The
  suffix must match exactly.
- Column names are whitespace-stripped on load to guard against stray spaces.
- Non-numeric E-values are coerced to NaN and dropped silently.
- The `1e-5` threshold is a loose presence call, it is not the KOfam GA
  threshold. KOfam threshold filtering is a separate upstream step.
- The output directory (`ATB_KO_Frequencies/`) must already exist or the write
  will fail with `FileNotFoundError`.

---

## `create_ko_counts.py`

Loads `ko_matrix.tsv`, re-binarises defensively (`> 0`), sums across columns
to get the number of genomes each KO appears in. Sorts descending.

**Output:** `ko_counts.tsv`: columns `KOs`, `KO_count`.

- Output index is named `KOs` (plural). `KO_frequency.py` looks for this name —
  don't rename it without updating that script.

---

## `create_ko_intersection.py`

Computes `M = X · Xᵀ` via sparse matrix multiplication. Entry (i, j) = number
of genomes containing both KO i and KO j. Diagonal M(i,i) = genome count for
KO i alone (equals `ko_counts`).

**Output:** `ko_intersection.tsv`: square symmetric KO × KO matrix.

**Nuances:**
- **Largest file in the pipeline** — at ~21,000 KOs the dense TSV output can
  be several GB. The sparse intermediate keeps memory manageable but the final
  write is dense.
- Does not binarise before multiplying (unlike `create_ko_counts.py`). If
  `ko_matrix.tsv` ever contains values > 1, intersection counts will be
  inflated. Worth adding `df = (df > 0).astype(np.uint8)` before the sparse
  cast for safety.
- The diagonal is the denominator for conditional probabilities in
  `build_neighbor_conditionals.py`: P(KO_j | KO_i) = M(i,j) / M(i,i).

---

## `run_ko_matrix.sh`

SLURM wrapper for `create_ko_matrix.py` only. Resources: 16 GB RAM, 4 CPUs
(unused script is single-threaded), 4-hour wall time, `twheeler` account,
`standard` partition.

---

**Fragile:** the conda init path (`/home/u13/nsontakke/miniconda3/...`) is
hardcoded to a specific user home directory.

## Portability

All paths in every script are hardcoded. If data moves, update these:

## Downstream

`ko_matrix.tsv` → `lineage_specific_ko_matrix.py` (splits by taxonomy)
`ko_intersection.tsv` → `build_neighbor_conditionals.py` (co-occurrence priors)
