# BLIMMP_2

# BLIMMP Removal Study — KO Recovery Benchmark

Benchmarks BLIMMP and six competing metabolic annotation tools on progressively
gene-depleted genomes to evaluate KO recovery as genome completeness decreases.

**Design:** 4 organisms × 11 removal levels (0–100%, step 10%) × 3 replicates = 132 samples.  
**Organisms:** *Acinetobacter baumannii*, *Prochlorococcus marinus* MED4, *P. marinus* MIT9313,
*Pseudomonas fluorescens* SBW25.  
**Tools compared:** HMMsearch (best-hit, E < 1e-5), KofamScan, anvi'o, METABOLIC, DRAM,
KEMET, BLIMMP (direct / neighbor / +sub, threshold 0.5).

All scripts run from `BASE_DIR = /xdisk/twheeler/nsontakke/Removal_Study_BLIMMP/` on puma.
Submit SLURM scripts with `sbatch` from that directory.

---

## Running order

```
1. run_tantan_prodigal.sh       # mask + call ORFs
2. run_busco.sh                 # QC (optional, independent)
3. run_subsample_orfs.sh        # create 11 × 3 ORF removal replicates
4. run_splice_orfs.sh           # mask removed ORF positions in genome
   └─ splice_orfs_from_genome.py
5. run_hmmsearch_spliced.sh     # HMMsearch on spliced .faa
6. run_kofamscan_spliced.sh     # KofamScan on spliced .fna
7. run_anvio_spliced.sh         # anvi'o on spliced .fna
8. run_metabolic.sh             # METABOLIC on spliced .fna
9. run_dram_spliced.sh          # DRAM on spliced .fna (filters masked contigs first)
   └─ filter_masked_contigs.py
10. run_kemet_spliced.sh        # KEMET on spliced .fna (needs KofamScan output)
11. fix_kofam_for_metapathpredict.py  # reformat KofamScan headers
12. run_new_matapatpredict.sh   # MetaPathPredict on fixed KofamScan files
13. parse_ko_presence_matrix.py # parse all tool outputs → binary KO matrices
14. parse_module_presence.py    # parse module-level presence (anvi'o, METABOLIC, DRAM)
15. parse_step_completeness.py  # parse step completeness (anvi'o, METABOLIC)
16. plot_*.py                   # generate figures (need KO_MATRICES/ ready)
```

---

## Scripts

### Genome preparation

**`run_tantan_prodigal.sh`**  
Masks repetitive regions with tantan, then runs prodigal to predict ORFs.  
Input: `{organism}/*.fna` — Output: `{organism}/{base}_masked.fna`, `{base}_ORFs.faa`  
`sbatch run_tantan_prodigal.sh`

**`run_busco.sh`**  
Runs BUSCO completeness assessment on reference genomes (QC step, not required downstream).  
Input: `{organism}/*.fna` — Output: BUSCO summary files per genome.  
`sbatch run_busco.sh`

**`subsample_orfs.py`**  
Randomly removes X% of ORFs from a `.faa` file across multiple percentages and replicates.  
Input: `-i ORFs.faa -n genome_name -o output_dir -p 0 10 20 ... 100 -r 3`  
Output: `{genome}_{pct}percentremoved_replicate{N}.faa` + `_removed.csv` per condition.

**`run_subsample_orfs.sh`**  
SLURM wrapper that calls `subsample_orfs.py` for each of the 4 reference genomes.  
Requires `run_tantan_prodigal.sh` to have completed first.  
`sbatch run_subsample_orfs.sh`

**`splice_orfs_from_genome.py`**  
Masks genomic positions of removed ORFs with N, producing a spliced genome `.fna`.  
Input: `-f genome.fna -a ORFs.faa -s subsample_dir -n genome_name -o output_dir`  
Output: `{sample}_spliced.fna` for each of the 132 samples.

**`run_splice_orfs.sh`**  
SLURM wrapper that calls `splice_orfs_from_genome.py` for all 4 genomes.  
Requires `run_subsample_orfs.sh` to have completed first.  
`sbatch run_splice_orfs.sh`

**`filter_masked_contigs.py`**  
Drops near-fully-N contigs from a spliced genome (DRAM's rRNA step fails on all-N contigs).  
Input: `-i spliced.fna -o filtered.fna -r report.tsv [-t 0.99]`  
Output: filtered `.fna` with all-N contigs removed. Called automatically by `run_dram_spliced.sh`.

---

### Tool execution

**`run_hmmsearch_spliced.sh`**  
Runs hmmsearch against 50 KO HMM chunks per sample, merges results into one `.domtblout`.  
Input: spliced `.faa` from `KOFAM_SPLICED/`  
Output: `HMMSEARCH_SPLICED/{sample}/{sample}_spliced.domtblout`  
`sbatch run_hmmsearch_spliced.sh`

**`run_kofamscan_spliced.sh`**  
Runs prodigal on each spliced `.fna` then annotates proteins with KofamScan (`--cut_tc`).  
Input: spliced `.fna` from `*_SPLICED_GENOMES/`  
Output: `KOFAM_SPLICED/{sample}/{sample}_kofam.tsv`  
`sbatch run_kofamscan_spliced.sh`

**`run_anvio_spliced.sh`**  
Builds anvi'o contigs DB, runs KOfam annotation, and estimates metabolism.  
Input: spliced `.fna` from `*_SPLICED_GENOMES/`  
Output: `ANVIO_SPLICED/{sample}/{sample}_hits.txt`, `_modules.txt`, `_metabolism.txt`  
`sbatch run_anvio_spliced.sh`

**`run_metabolic.sh`**  
Runs METABOLIC-G on each spliced genome.  
Input: spliced `.fna` from `*_SPLICED_GENOMES/`  
Output: `METABOLIC_SPLICED/METABOLIC_{sample}/KEGG_identifier_result/total.hits.txt`  
`sbatch run_metabolic.sh`

**`run_dram_spliced.sh`**  
Filters masked contigs with `filter_masked_contigs.py`, then runs DRAM annotate + distill.  
Input: spliced `.fna` from `*_SPLICED_GENOMES/`  
Output: `DRAM_SPLICED/{sample}/annotations.tsv`  
`sbatch run_dram_spliced.sh`

**`run_kemet_spliced.sh`**  
Runs KEMET using KofamScan annotations as input.  
Input: spliced `.fna` + `KOFAM_SPLICED/{sample}/{sample}_kofam.tsv`  
Output: `KEMET_SPLICED/{sample}/ktests/*.ktest`  
Requires `run_kofamscan_spliced.sh` to have completed first.  
`sbatch run_kemet_spliced.sh`

**`fix_kofam_for_metapathpredict.py`**  
Reformats KofamScan `.tsv` column headers to match MetaPathPredict's expected format.  
Input: `KOFAM_SPLICED/` — Output: `KOFAM_SPLICED_FIXED/` (originals unchanged)  
`conda activate blimmp-work && python scripts/fix_kofam_for_metapathpredict.py`

**`run_new_matapatpredict.sh`**  
Runs MetaPathPredict on all fixed KofamScan files to predict module presence.  
Input: `KOFAM_SPLICED_FIXED/**/*_kofam.tsv`  
Output: `METAPATHPREDICT_SPLICED/metapathpredict_all.tsv`  
Requires `fix_kofam_for_metapathpredict.py` first.  
`sbatch run_new_matapatpredict.sh`

**`remove_metab.sh`**  
Utility: deletes `METABOLIC_SUBSAMPLES/` to free disk space.  
`sbatch remove_metab.sh`

---

### Parsing (run after all tools complete)

**`parse_ko_presence_matrix.py`**  
Parses all 9 tool outputs into binary KO presence/absence matrices restricted to the
bacterial KEGG module KO universe (`module_ko_reaction.json`).  
Input: all tool output directories under `BASE_DIR/`  
Output: `KO_MATRICES/sample_metadata.csv` + one `*_binary.csv` per tool  
`conda activate blimmp-work && python scripts/parse_ko_presence_matrix.py`

**`parse_module_presence.py`**  
Aggregates module-level presence calls from anvi'o, METABOLIC, and DRAM.  
Input: `ANVIO_SPLICED/`, `METABOLIC_SPLICED/`, `DRAM_SPLICED/`  
Output: `module_presence.csv`  
`conda activate blimmp-work && python scripts/parse_module_presence.py`

**`parse_step_completeness.py`**  
Extracts per-step module completeness scores from anvi'o and METABOLIC.  
Input: `ANVIO_SPLICED/` (`*_modules.txt`) and `METABOLIC_SPLICED/` (`worksheet4.tsv`)  
Output: `step_completeness.csv`  
`conda activate blimmp-work && python scripts/parse_step_completeness.py`

---

### Visualization (run after `parse_ko_presence_matrix.py`)

All plot scripts read from `KO_MATRICES/` and write to `PLOTS/`.  
`conda activate blimmp-work && python scripts/<script>.py`

**`plot_ko_recovery_100pct.py`**  
Self-recovery: fraction of each tool's own 0%-removal KOs still detected at higher removal levels.  
Output: `PLOTS/ko_recovery_100pct.{png,pdf}`

**`plot_ko_kofamscan_recovery.py`**  
Recall and precision vs KofamScan stable KOs at 0% removal. Layout: 2 rows × 4 organism columns.  
Output: `PLOTS/ko_kofamscan_recovery.{png,pdf}`

**`plot_ko_hmmsearch_recovery.py`**  
Recall and precision vs HMMsearch (E < 1e-5, best hit per ORF) stable KOs at 0% removal.
X-axis truncated at 90% (100% removal is trivially empty).  
Output: `PLOTS/ko_hmmsearch_recovery.{png,pdf}`

**`plot_ko_quadrant.py`**  
Scatter plot: BLIMMP `ko_probability` (X) vs. number of non-BLIMMP tools detecting each KO (Y),
at 0% removal. Four quadrants at X = 0.5, Y = 0.5.  
Reads raw BLIMMP output files in addition to `KO_MATRICES/`.  
Output: `PLOTS/ko_quadrant.{png,pdf}`

---

## Environment

```bash
conda activate blimmp-work   # Python ≥ 3.9, pandas, numpy, matplotlib
```
