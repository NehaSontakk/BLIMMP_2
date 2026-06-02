# HMMER Adjudication

Scripts for processing HMMER `.domtblout` output files into one clean KO annotation per genomic locus, for use in building reference genome priors.

Both this pipeline and BLIMMP use identical adjudication logic to decide which KO wins each genomic locus:

1. Parse `.domtblout`, auto-detect KO column
2. Compute HMM union coverage per (strand, target, KO) group via Numba
3. Cluster overlapping ORF hits by genomic position (60% fractional overlap threshold)
4. Compute softmax within each overlap group
5. Apply noise-floor confidence (E-value threshold = 1e-4) to pick one winner per overlap group
6. Genome-wide KO dedup — keep top-scoring hit per KO

### The adjudication process differs from BLIMMP pipeline in the use of KOFAM threshold.

| | This pipeline | BLIMMP |
|---|---|---|
| Goal | Hard winner per locus for reference genomes | Soft probability per KO for query genome inference |
| KOfam threshold | Not applied | Applied — recalculates hit_conf |
| KO universe | All hmmsearch hits | Module-relevant KOs only |
| Output | One row per unique KO | Full scored DataFrame with flags |

The adjudication scores are identical for matched KOs. See `Comparison_with_BLIMMP/` for validation.

---

## Usage

### Step 1 — Generate commands

```bash
python commands_preprocess.py
```

Scans `ATB_Dechunked_HMMER/` for `.domtblout` files and writes one preprocess command per genome to `commands.txt`. Core adjudication script, processes one `.domtblout` into one winner per KO.

### Step 2 — Run adjudication (parallel)

```bash
sbatch parallel_preprocess.sh
```

Runs 30 genomes per SLURM task in parallel using GNU parallel. Adjust `--array` range based on total number of genomes (total commands / 30, rounded up).

### Step 2 (alternative) — Run serially

```bash
sbatch run_preprocess_domtblout.sh
```

Generates `commands.txt` one preprocess command per genome created as a fallback. Processes all genomes sequentially. Skips already-completed files so safe to rerun for recovery.

---

## Output

One CSV per genome in `PROCESS_HMMER_OUTPUT/`:

```
{sample}_grouped_hits_dedup.csv
```

Each row is one unique KO detected in the genome, with columns:

| Column | Description |
|---|---|
| `KO id` | KEGG Orthology identifier |
| `target name` | ORF name |
| `score` | Full-sequence HMM score |
| `hit_conf` | Noise-floor annotation confidence (0–1) |
| `overlap_group` | Genomic locus identifier |
| `hmm_coverage_fraction` | Fraction of HMM covered by alignment |

---

## Hardcoded Paths — Update Before Running

| File | Path to update |
|---|---|
| `commands_preprocess.py` | `INDIR` — path to `.domtblout` files |
| `commands_preprocess.py` | `OUTDIR` — path for output CSVs |
| `run_preprocess_domtblout.sh` | `INDIR` and `OUTDIR` |
| Both `.sh` scripts | `--account=twheeler` → your SLURM account |
| Both `.sh` scripts | `conda activate test` → your conda environment |

---

## Notes

- `parallel_preprocess.sh` is a sample script — the original run used array range `663-714` as a recovery pass. For a fresh run use `--array=1-N` where N = ceil(total_genomes / 30).
- Requires `numba`, `pandas`, `numpy` — install via conda or pip.
