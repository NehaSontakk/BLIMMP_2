# 1. Tantan and Prodigal Scripts

Scripts for repeat-masking genome assemblies and predicting open reading frames (ORFs).


## Step 1 — Tantan (Repeat Masking)

### Script
`Tantan_Scripts/run_tantan.sh`


Runs [Tantan](https://gitlab.com/mcfrith/tantan) on each raw genome `.fa` file to mask repetitive regions (replacing them with `N`s). This reduces spurious ORF predictions in repetitive regions downstream.

### Input
`ATB_FASTA/*.fa` — raw genome assemblies (produced by unpacking AllTheBacteria downloads)

### Output
`ATB_TANTAN/*_masked.fa` — repeat-masked genome assemblies

### Usage
```bash
sbatch Tantan_Scripts/run_tantan.sh
```

### Notes
- Runs as a SLURM array job (44 tasks × 500 files per task)
- Skips files that are already masked
- Requires Tantan installed, update this path for your environment.

---

## Step 2 — Prodigal (ORF Prediction)

### Scripts
- `Prodigal_Scripts/generate_prodigal_commands.py`
- `Prodigal_Scripts/run_prodigal.sh`

Predicts protein-coding ORFs in each masked genome using [Prodigal](https://github.com/hyattpd/Prodigal).

### Input
`ATB_TANTAN/*_masked.fa` repeat-masked genome assemblies from Step 1

### Output
`ATB_PRODIGAL/` — for each genome:
- `*_ORFs.faa` — predicted protein sequences
- `*_ORFs.fna` — predicted nucleotide sequences
- `*_ORFs.gff` — gene coordinates

### Usage

**1. Generate commands:**
```bash
python Prodigal_Scripts/generate_prodigal_commands.py
```
This scans `ATB_TANTAN/` for masked files, skips any already processed, and writes `prodigal_cmds.txt`.

**2. Run Prodigal:**
```bash
sbatch Prodigal_Scripts/run_prodigal.sh
```

### Notes
- `generate_prodigal_commands.py` is incremental — safe to rerun if jobs fail; it skips already-completed samples
- SLURM array runs up to 500 tasks, each processing a chunk of commands from `prodigal_cmds.txt`
- Requires Prodigal available on your `$PATH` or update the path in the script


# 2. HMMER Scripts

Scripts for running hmmsearch against the KEGG KO HMM database across all genome ORFs, detecting and recovering missing chunks, and combining outputs.

## Overview

```
Use the Prodigal *_ORFs.faa
    → step1_detect_missing_chunks.py  (find missing chunk/sample pairs)
    → step2_hmmsearch_missing.py      (generate hmmsearch commands)
COMMANDS/commands_1..130.txt
    → run_hmmsearch.sh                (run hmmsearch in parallel)
ATB_HMMER_Run2/{sample}_chunk{N}.domtblout
ATB_HMMER_Run2/{sample}_chunk{N}.hmmout
    → combine_domtblout.py            (merge chunked .domtblout per sample)
    → combine_hmmout.py / combine_hmmout.sh  (merge chunked .hmmout per sample)
ATB_Dechunked_HMMER/{sample}.combined.domtblout   final output
ATB_Dechunked_HMMER_hmmout/{sample}.hmmout        final output
```

## Database

The HMM database is the **KOfam (KEGG Orthology HMM)** database, split into 130 chunks
(`chunked_hmmdb/chunk_1.hmm` ... `chunk_130.hmm`). Each genome is searched against all
130 chunks, producing 130 output files per genome.

Download KOfam:
```bash
wget ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz
wget ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz
tar xf profiles.tar.gz
gunzip ko_list.gz
```

The profiles then need to be split into 130 chunks for parallel searching. See
`step2_hmmsearch_missing.py` for how chunk indices map to HMM files.

---

## Usage

### Step 1 — Detect missing chunks

```bash
python step1_detect_missing_chunks.py
```

Scans `ATB_HMMER_Run2/` for completed `.domtblout` files and compares against the expected 130 chunks × all genomes. On the first run with an empty output directory, everything will be flagged as missing — this is expected. Writes `missing_chunks.txt`.

### Step 2 — Generate hmmsearch commands

```bash
python step2_hmmsearch_missing.py
```

Reads `missing_chunks.txt` and writes one hmmsearch command per missing chunk to `COMMANDS/commands_1..130.txt`. Each command searches one genome's ORFs against one HMM chunk.

### Step 3 — Run hmmsearch

```bash
sbatch run_hmmsearch.sh
```

Runs hmmsearch in parallel across all 130 chunks as a SLURM array job. Each task processes one `COMMANDS/commands_N.txt` file using GNU parallel.

To rerun only a subset of chunks:
```bash
sbatch --array=1-50 run_hmmsearch.sh
```

### Step 4 — Combine outputs

```bash
# Combine .hmmout files
sbatch combine_hmmout.sh

# Combine .domtblout files
python combine_domtblout.py -i ATB_HMMER_Run2/
```

---

## Hardcoded Paths need update before running

| File | Hardcoded Path | What to update |
|---|---|---|
| `step1_detect_missing_chunks.py` | `/xdisk/cgoubert/nsontakke/ATB_HMMER_Run2/` | Your hmmsearch output directory |
| `step1_detect_missing_chunks.py` | `/xdisk/twheeler/nsontakke/ATB_Analysis_0725/ATB_PRODIGAL/` | Your Prodigal output directory |
| `step2_hmmsearch_missing.py` | `/xdisk/cgoubert/nsontakke/ATB_HMMER_Run2` | Your hmmsearch output directory |
| `step2_hmmsearch_missing.py` | `/xdisk/twheeler/nsontakke/ATB_Analysis_0725/ATB_PRODIGAL` | Your Prodigal output directory |
| `step2_hmmsearch_missing.py` | `/xdisk/twheeler/nsontakke/ATB_Analysis_0725/chunked_hmmdb` | Your KEGG HMM database directory |
| `combine_hmmout.py` | `ATB_HMMER_Run2` | Must match your hmmsearch output directory |
| `run_hmmsearch.sh` | `OUTPUT_DIR="/path/to/your/ATB_HMMER_Run2"` | Your hmmsearch output directory |
| `run_hmmsearch.sh` | `--account=twheeler` | Your SLURM account |

---

## Notes

- `run_hmmsearch.sh` is a **sample script** — the original runs were split across two clusters (elgato, puma) and two accounts (twheeler, cgoubert) due to resource constraints. This consolidated script covers all 130 chunks in one job.
- Each genome produces 130 chunk output files before combining — ensure sufficient disk space.
- `combine_domtblout.py` expects chunked files named `{sample}_chunk{N}.domtblout` in the input directory.
