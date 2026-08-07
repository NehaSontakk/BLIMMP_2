# Removal_Study_BLIMMP

ORF-removal ablation study: for each of several bacterial/cyanobacterial genomes,
progressively remove a percentage of predicted ORFs and measure how well four
different metabolic-annotation tools (BLIMMP, anvi'o, DRAM, METABOLIC) recover
their own full-genome KO/module/pathway calls as information is degraded.

## Pipeline overview

```
tantan (repeat masking) + prodigal (gene calling)
        |
        +-- BUSCO (genome completeness, independent check)
        |
        +-- HMMER: split all_profiles.hmm -> parallel hmmsearch -> merge
        |          (full-genome KO search, ground truth for BLIMMP)
        |
        +-- ORF subsampling (0/10/20/.../100% removed x 3 replicates)
                |
                +-- DRAM: annotate_genes -> distill
                +-- anvi'o: gen-contigs-database -> run-kegg-kofams -> estimate-metabolism
                +-- METABOLIC: METABOLIC-G.pl
```

Genomes currently included: `Acinetobacter_baumannii`, `MED4`, `MIT9313`, `SS120`,
`Pseudomonas_fluorescens_SBW25`.

## Directory layout

```
<Genome>/<Genome>_SUBSAMPLES/          per-genome subsampled ORF .faa + removed-ID .csv files
DRAM_SUBSAMPLES/<sample>/              DRAM annotate_genes + distill output, per sample
ANVIO_SUBSAMPLES/<sample>/             anvi'o contigs.db + estimate-metabolism output, per sample
METABOLIC_SUBSAMPLES/METABOLIC_<sample>/   METABOLIC-G.pl output, per sample
HMMER_FULL_GENOME/<genome>/            merged full-genome HMMER domtblout/tblout
HMM_CHUNKS/                            all_profiles.hmm split into 50 chunks for parallel search
logs/                                  all SLURM job logs (*_%A_%a.out/.err)
.make_stamps/                          Makefile bookkeeping -- do not edit, see below
```

## Running the pipeline

Everything is orchestrated by the `Makefile`, which submits each stage as a SLURM
job via `sbatch` and blocks until it finishes before moving to the next stage.
Each stage only reruns if it hasn't completed before (tracked via stamp files in
`.make_stamps/`), so `make` is safe to rerun after an interruption.

```bash
make            # run everything not yet done, in dependency order
make status     # show which stages are complete
make subsample  # run only up through subsampling
make metabolic  # run only up through METABOLIC (and its prerequisites)
make clean-hmmer   # force the HMMER stage (and re-dependents) to rerun
make clean         # wipe all progress tracking, start over
```

**`make` blocks and polls until each SLURM job finishes**, which can take hours
across the full pipeline. Run it inside `tmux`/`screen`, or as its own lightweight
background job, so it survives a disconnect:

```bash
tmux new -s pipeline
make
# Ctrl-b d to detach; tmux attach -t pipeline to check back in
```

See comments at the top of `Makefile` for the full stage list and dependency graph.

## Scripts

| Script | Purpose |
|---|---|
| `run_tantan_prodigal.sh` | Repeat-mask genomes (tantan), call ORFs (prodigal) |
| `run_busco.sh` | BUSCO completeness per genome |
| `subsample_orfs.py` / `run_subsample_orfs.sh` | Random ORF removal at 0-100%, 3 replicates, with built-in consistency checks |
| `split_hmm_profiles.sh` / `run_hmmsearch_array.sh` / `merge_hmmsearch_results.sh` | Full-genome HMMER search, chunked for parallelism |
| `run_dram_array.sh` / `run_dram_distill_array.sh` | DRAM annotation + distillation per subsample |
| `make_external_gene_calls.py` / `run_anvio_build_array.sh` / `run_anvio_estimate_array.sh` | anvi'o contigs-db build + KO/module/path estimation per subsample |
| `run_metabolic_array.sh` | METABOLIC-G.pl per subsample |
| `finish_anvio.sh` | One-time KEGG data download for anvi'o (`anvi-setup-kegg-data`) |
| `Makefile` | Orchestrates all of the above in dependency order |

## Requirements

- SLURM cluster access (`--account=twheeler`, `--partition=standard`)
- Conda envs: `my_dram_env`, `anvio-9`, `METABOLIC_v4.0`
- KEGG data set up once via `finish_anvio.sh` before running anvi'o stages
