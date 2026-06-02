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
