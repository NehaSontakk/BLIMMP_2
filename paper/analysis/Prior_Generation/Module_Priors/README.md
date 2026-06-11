# Module_Priors

Generates per-module presence frequencies from METABOLIC outputs across all ATB genomes.
The final output (`module_freq.txt`) is used directly by BLIMMP as the Beta prior mean for each module.

## Scripts

| Script | Input | Output | Description |
|---|---|---|---|
| `generate_metabolic_commands.py` | `ATB_PRODIGAL/*.faa` | `metabolic_commands.txt` | Creates one METABOLIC-G command per genome; copies `.faa` into per-genome subdirectory under `ATB_METABOLIC/` |
| `metabolic_batch_creation.sh` | `metabolic_commands.txt` | `metabolic_batches/batch_1.txt` … `batch_100.txt` | Splits 21,392 commands evenly across 100 batch files for SLURM array submission |
| `metabolic_run_array.sh` | `metabolic_batches/batch_${SLURM_ARRAY_TASK_ID}.txt` | `ATB_METABOLIC/*_ORFs_METABOLIC/` | SLURM array (1–100); runs each batch sequentially, logs per-command success/failure |
| `metabolic_files.sh` | `ATB_METABOLIC/*/METABOLIC_result_each_spreadsheet/METABOLIC_result_worksheet3.tsv` | `metabolic_sample_module_lastcol.tsv` | Scrapes module ID and presence call from worksheet3 of every genome; outputs one row per (genome, module) |
| `count_modules.py` | `metabolic_sample_module_lastcol.tsv` | `module_freq.txt`, `module_presence_fraction.tsv` | Counts Present/Absent per module across all genomes; divides by total unique genomes to get frequency |

---

## Inputs

- `ATB_PRODIGAL/*.faa` — Prodigal-predicted ORFs, one `.faa` per genome (21,392 files)
- METABOLIC v4.0 installed at `/xdisk/twheeler/nsontakke/Software/METABOLIC_running_folder/METABOLIC/`

## Outputs

- `metabolic_sample_module_lastcol.tsv` — combined presence/absence table (Sample_ID, module_id, presence); ~10M rows; not versioned in git
- `module_freq.txt` — tab-separated (module_id, frequency); direct input to BLIMMP prior parameterisation
- `module_presence_fraction.tsv` — same frequencies with raw present counts included

---

## Nuances

worksheet3 is the correct file.

**Denominator is total unique genomes, not row count.** `count_modules.py` divides by `n_genomes = df["Sample_ID"].nunique()` (21,391), not by the per-module row count. These differ if any genome is missing a module entry, which would inflate frequency estimates.

**`metabolic_commands.txt` is ephemeral.** It is overwritten each time `generate_metabolic_commands.py` runs. The batch files in `metabolic_batches/` are the durable record of what ran.
