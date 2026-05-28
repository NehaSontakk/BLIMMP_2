# ATB Download Scripts

Scripts for downloading and unpacking bacterial genome assemblies from [AllTheBacteria](https://allthebacteria.org/).

## Overview

This pipeline filters high-quality bacterial assemblies from AllTheBacteria, downloads them in batches via SLURM, and unpacks them for downstream analysis. It selects up to 5 samples per species, resulting in ~21,878 genomes.

## Files

| File | Description |
|---|---|
| `sample_and_prep_ATB_downloads.py` | Filters samples and generates download/unzip command files |
| `filtered_samples_filelist_080725.tsv` | The exact sample list used in this study (August 7, 2025) |
| `download_files.sh` | SLURM array job to download `.tar.xz` archives from OSF |
| `unpack_downloads.sh` | SLURM array job to unpack the downloaded archives |

## Prerequisites

### 1. Download the AllTheBacteria file list

```bash
wget https://osf.io/download/4yv85/ -O file_list.all.latest.tsv.gz
gunzip file_list.all.latest.tsv.gz
```

> **Note:** The URL above may change. If it fails, check the current URL at https://osf.io/zxfmy/files/osfstorage and look for `file_list.all.latest.tsv.gz`.

### 2. Download the high-quality sample list

```bash
wget https://osf.io/download/hq_set.sample_list.txt.gz -O hq_set.sample_list.txt.gz
gunzip hq_set.sample_list.txt.gz
mv hq_set.sample_list.txt hq_seq.samples.txt
```

> If the link above fails, download `hq_set.sample_list.txt.gz` manually from https://osf.io/zxfmy/files/osfstorage.

## Usage

### Step 1 — Filter samples and generate download commands

```bash
python sample_and_prep_ATB_downloads.py
```

This reads `hq_seq.samples.txt` and `file_list.all.latest.tsv`, filters to up to 5 samples per species, and outputs:
- `filtered_samples_filelist_080725.tsv` — the filtered sample list
- `download_commands.txt` — wget commands for each archive
- `unzip_commands.txt` — tar commands for each archive

> To reproduce the exact samples used in this study, skip this step and use the `filtered_samples_filelist_080725.tsv` already provided in this folder.

### Step 2 — Download archives (SLURM)

```bash
sbatch download_files.sh
```

Downloads `.tar.xz` archives into `ATB_Downloads/`. Skips files that already exist, so it's safe to rerun if interrupted.

### Step 3 — Unpack archives (SLURM)

```bash
sbatch unpack_downloads.sh
```

Unpacks each archive in batches of 500.

## Output

Unpacked `.fa` genome files in subdirectories under `ATB_Downloads/`.

## Notes

- Logs are written to `logs/` which is created automatically.
- `download_commands.txt` and `unzip_commands.txt` are generated files and are not stored in this repo.
