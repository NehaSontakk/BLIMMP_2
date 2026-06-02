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
