#!/bin/bash
# Run this ON PUMA from inside /xdisk/twheeler/nsontakke/Removal_Study_BLIMMP
#
# Creates an ORPHAN branch on BLIMMP_2 -- i.e. it does NOT start from main's
# history or files. The branch will contain ONLY what's copied in below:
# the code files and small result files from this directory. Data/output
# subdirectories (genomes, subsample results, HMM chunks, logs) are skipped.

set -euo pipefail

REPO_DIR="/xdisk/twheeler/nsontakke/Removal_Study_BLIMMP/BLIMMP_2"  # where BLIMMP_2 is cloned
BRANCH_NAME="removal-study-orf-subsampling"     # change if you want a different name
SRC_DIR="/xdisk/twheeler/nsontakke/Removal_Study_BLIMMP"            # the directory with your scripts

# 1. Clone the repo if it isn't already on Puma
if [ ! -d "$REPO_DIR" ]; then
  echo "Cloning BLIMMP_2 into $REPO_DIR ..."
  git clone git@github.com:NehaSontakk/BLIMMP_2.git "$REPO_DIR"
fi

cd "$REPO_DIR"
git fetch origin

# 2. Create an orphan branch (no history, no files from main)
if git show-ref --verify --quiet "refs/heads/$BRANCH_NAME"; then
  echo "Local branch $BRANCH_NAME already exists -- checking it out."
  git checkout "$BRANCH_NAME"
else
  git checkout --orphan "$BRANCH_NAME"
fi

# 3. Clear out anything currently tracked/staged (orphan branch starts with
#    main's working files still on disk -- we don't want those)
git rm -rf --cached . >/dev/null 2>&1 || true
find . -mindepth 1 -not -path './.git*' -delete

# 4. Copy over just the code + small result files (not data/output dirs)
cp "$SRC_DIR"/*.py .  2>/dev/null || true
cp "$SRC_DIR"/*.sh .  2>/dev/null || true
cp "$SRC_DIR"/completeness_stats.csv . 2>/dev/null || true
cp "$SRC_DIR"/Makefile . 2>/dev/null || true
cp "$SRC_DIR"/README* . 2>/dev/null || true

# NOTE: deliberately NOT copying these (data/output directories):
#   Acinetobacter_baumannii, ANVIO_SUBSAMPLES, DRAM_SUBSAMPLES, HMM_CHUNKS,
#   MED4, METABOLIC_SUBSAMPLES, MIT9313, Pseudomonas_fluorescens_SBW25,
#   SS120, logs

# 5. Commit and push
git add -A
git status   # sanity check before committing -- review what's staged
git commit -m "ORF subsampling removal study: code + summary stats (Acinetobacter, Pseudomonas SBW25, Prochlorococcus ecotypes)"
git push -u origin "$BRANCH_NAME"

echo ""
echo "Done. Orphan branch '$BRANCH_NAME' pushed to origin."
echo "View it at: https://github.com/NehaSontakk/BLIMMP_2/tree/$BRANCH_NAME"
