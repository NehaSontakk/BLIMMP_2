# Tantan and Prodigal Scripts

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
