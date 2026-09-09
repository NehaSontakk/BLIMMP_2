#!/usr/bin/env python3
"""
Step 1 of 2 — Parse all tool outputs into binary KO presence/absence matrices.

Run ONCE on puma (slow — reads thousands of files).
Outputs go to KO_MATRICES/ under BASE_DIR:

  KO_MATRICES/
    sample_metadata.csv          — sample x (organism, removal_pct, replicate)
    HMMsearch_binary.csv         — sample x KO (0/1), raw HMM hits E-value < 1e-5
    KofamScan_binary.csv         — KofamScan star-filter passing hits
    anvio_binary.csv
    METABOLIC_binary.csv
    DRAM_binary.csv
    KEMET_binary.csv
    BLIMMP_prior_binary.csv      — pass-1 prior: direct HMM evidence, ko_probability_prior > 0.5
    BLIMMP_nosub_binary.csv      — pass-1 posterior: direct + neighbor context, ko_probability > 0.5
    BLIMMP_binary.csv            — pass-2: direct + neighbor + substitute, ko_probability > 0.5

IMPORTANT: All matrices are restricted to the bacterial KEGG module KO universe
(KOs appearing in BLIMMP_8Sep2026/BLIMMP/module_ko_reaction.json). This ensures
a fair comparison — BLIMMP only searches these KOs, so all tools are evaluated
on the same set.

Usage:
    conda activate blimmp-work
    python parse_ko_matrices.py

MetaPathPredict is excluded — it outputs module predictions, not KO lists.
"""

import json
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE_DIR  = Path("/xdisk/twheeler/nsontakke/Removal_Study_BLIMMP")
OUT_DIR   = BASE_DIR / "KO_MATRICES"
OUT_DIR.mkdir(exist_ok=True)

# ---------------------------------------------------------------------------
# Bacterial KEGG module KO universe
# BLIMMP only searches KOs in bacterial modules; restrict all tools to this
# universe so the precision/recall comparison is apples-to-apples.
# ---------------------------------------------------------------------------
MODULE_KO_JSON = BASE_DIR / "BLIMMP_8Sep2026" / "BLIMMP" / "module_ko_reaction.json"
with open(MODULE_KO_JSON) as _f:
    _module_data = json.load(_f)

BACTERIAL_MODULE_KOS: set = set()
for _module_id, _ko_dict in _module_data.items():
    BACTERIAL_MODULE_KOS.update(_ko_dict.keys())

print(f"Bacterial KEGG module KOs loaded from module_ko_reaction.json: "
      f"{len(BACTERIAL_MODULE_KOS):,} unique KOs")
print()

# ---------------------------------------------------------------------------
# Sample name parser
# Pattern: {organism}_{pct}percentremoved_replicate{N}
# ---------------------------------------------------------------------------
SAMPLE_RE = re.compile(r"^(.+)_(\d+)percentremoved_replicate(\d+)$")

def parse_sample_name(name):
    m = SAMPLE_RE.match(name)
    return (m.group(1), int(m.group(2)), int(m.group(3))) if m else None

# ---------------------------------------------------------------------------
# Per-tool KO set parsers
# Each returns a set of KO strings (e.g. {"K00001", "K00134", ...})
# ---------------------------------------------------------------------------

def parse_hmmsearch(sample: str) -> set:
    """*_spliced.domtblout: best-hit per ORF with full-sequence E-value < 1e-5.

    Each ORF (query protein) is assigned to the single KO HMM with the lowest
    E-value among all hits passing E < 1e-5.  This avoids counting multiple KOs
    from weak, overlapping hits on the same ORF.

    Handles both hmmscan (KO in col 0, ORF in col 2) and hmmsearch
    (KO in col 3, ORF in col 0) layouts.
    """
    p = BASE_DIR / "HMMSEARCH_SPLICED" / sample / f"{sample}_spliced.domtblout"
    if not p.exists():
        return set()
    ko_re = re.compile(r"^K\d{5}$")
    best_hits: dict = {}   # orf → (best_evalue, best_ko)
    with open(p) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 7:
                continue
            try:
                evalue = float(fields[6])   # full-sequence E-value (col 6)
            except ValueError:
                continue
            if evalue >= 1e-5:
                continue
            # Determine KO and ORF based on layout
            if ko_re.match(fields[0]):
                ko, orf = fields[0], fields[2]   # hmmscan
            elif len(fields) > 3 and ko_re.match(fields[3]):
                ko, orf = fields[3], fields[0]   # hmmsearch
            else:
                continue
            if orf not in best_hits or evalue < best_hits[orf][0]:
                best_hits[orf] = (evalue, ko)
    return {ko for _, ko in best_hits.values()}


def parse_blimmp_prior(sample: str) -> set:
    """*__BLIMMP_dk.csv: ko_probability_prior > 0.5 (pass-1 prior, direct HMM evidence only)."""
    p = BASE_DIR / "BLIMMP_SPLICED" / sample / f"{sample}__BLIMMP_dk.csv"
    if not p.exists():
        return set()
    df = pd.read_csv(p, usecols=["ko_id", "ko_probability_prior"])
    return set(df.loc[df["ko_probability_prior"] > 0.5, "ko_id"].str.strip())


def parse_blimmp_nosub(sample: str) -> set:
    """*__BLIMMP_dk.csv: ko_probability > 0.5 (pass-1 posterior, direct + neighbor context)."""
    p = BASE_DIR / "BLIMMP_SPLICED" / sample / f"{sample}__BLIMMP_dk.csv"
    if not p.exists():
        return set()
    df = pd.read_csv(p, usecols=["ko_id", "ko_probability"])
    return set(df.loc[df["ko_probability"] > 0.5, "ko_id"].str.strip())


def parse_blimmp(sample: str) -> set:
    """Full BLIMMP pipeline: pass-1 neighbor KOs ∪ pass-2 substituted KOs.

    Takes the union of parse_blimmp_nosub() and the substituted_dk KOs so
    that BLIMMP is always a superset of BLIMMP_nosub — substitutions can only
    add KOs, never remove them.  Both use threshold > 0.5.
    """
    base = parse_blimmp_nosub(sample)
    d = BASE_DIR / "BLIMMP_SPLICED" / sample
    hits = list(d.glob("*_substituted_dk.csv")) if d.exists() else []
    if not hits:
        return base
    df = pd.read_csv(hits[0], usecols=["ko_id", "ko_probability"])
    sub_kos = set(df.loc[df["ko_probability"] > 0.5, "ko_id"].str.strip())
    return base | sub_kos


def parse_kofamscan(sample: str) -> set:
    """*_kofam.tsv (tab-separated): lines where field[0] == '*' are passing hits.
    Fields: *  <gene>  <KO>  <threshold>  <score>  <evalue>  <definition>
    KO is at field index 2."""
    p = BASE_DIR / "KOFAM_SPLICED" / sample / f"{sample}_kofam.tsv"
    kos = set()
    if not p.exists():
        return kos
    with open(p, "rb") as f:
        for raw in f:
            row = raw.decode(errors="replace")
            fields = row.split("\t")
            if fields[0] == "*" and len(fields) >= 3:
                kos.add(fields[2].strip())
    return kos


def parse_anvio(sample: str) -> set:
    """*_hits.txt: 'enzyme' column. anvi'o already applies its score threshold."""
    p = BASE_DIR / "ANVIO_SPLICED" / sample / f"{sample}_hits.txt"
    if not p.exists():
        return set()
    df = pd.read_csv(p, sep="\t", usecols=["enzyme"])
    kos = df["enzyme"].dropna().astype(str).str.strip()
    return {k for k in kos if re.match(r"^K\d{5}$", k)}


def parse_metabolic(sample: str) -> set:
    """total.hits.txt: KO detected if line has >= 2 whitespace-separated tokens.
    'K00003 NZ_CP091367.1_3055_1' → detected; 'K00001' alone → absent."""
    p = (BASE_DIR / "METABOLIC_SPLICED"
         / f"METABOLIC_{sample}"
         / "KEGG_identifier_result"
         / "total.hits.txt")
    kos = set()
    if not p.exists():
        return kos
    with open(p) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                kos.add(parts[0].strip())
    return kos


def parse_dram(sample: str) -> set:
    """annotations.tsv: 'ko_id' column; non-empty K-number entries."""
    p = BASE_DIR / "DRAM_SPLICED" / sample / "annotations.tsv"
    if not p.exists():
        return set()
    df = pd.read_csv(p, sep="\t", usecols=["ko_id"])
    kos = df["ko_id"].dropna().astype(str).str.strip()
    return {k for k in kos if re.match(r"^K\d{5}$", k)}


def parse_kemet(sample: str) -> set:
    """ktests/*.ktest: one KO per line, already filtered by KEMET."""
    ktest_dir = BASE_DIR / "KEMET_SPLICED" / sample / "ktests"
    hits = list(ktest_dir.glob("*.ktest")) if ktest_dir.exists() else []
    if not hits:
        return set()
    kos = set()
    with open(hits[0]) as f:
        for line in f:
            k = line.strip()
            if k:
                kos.add(k)
    return kos


PARSERS = {
    "HMMsearch":     parse_hmmsearch,     # raw HMM hits, E-value < 1e-5
    "KofamScan":     parse_kofamscan,     # HMM with family-specific score threshold
    "anvio":         parse_anvio,         # filesystem-safe name (no apostrophe)
    "METABOLIC":     parse_metabolic,
    "DRAM":          parse_dram,
    "KEMET":         parse_kemet,
    "BLIMMP_prior":  parse_blimmp_prior,  # pass-1 prior: direct HMM evidence, ko_probability_prior > 0.01
    "BLIMMP_nosub":  parse_blimmp_nosub,  # pass-1 posterior: direct + neighbor context, ko_probability > 0.01
    "BLIMMP":        parse_blimmp,        # pass-2: direct + neighbor + substitute, ko_probability > 0.25
}

# ---------------------------------------------------------------------------
# Discover samples
# ---------------------------------------------------------------------------
print("Discovering samples from KOFAM_SPLICED ...")
kofam_root = BASE_DIR / "KOFAM_SPLICED"
all_dirs   = sorted(d.name for d in kofam_root.iterdir() if d.is_dir())

samples  = []
meta_rows = []
for name in all_dirs:
    parsed = parse_sample_name(name)
    if parsed is None:
        print(f"  SKIP (cannot parse): {name}")
        continue
    organism, pct, rep = parsed
    samples.append(name)
    meta_rows.append(dict(sample=name, organism=organism,
                          removal_pct=pct, replicate=rep))

print(f"  {len(samples)} samples")

# Save metadata
meta_df = pd.DataFrame(meta_rows).set_index("sample")
meta_df.to_csv(OUT_DIR / "sample_metadata.csv")
print(f"  Saved: {OUT_DIR}/sample_metadata.csv")
print()

# ---------------------------------------------------------------------------
# Parse each tool → build binary matrix → save CSV
# ---------------------------------------------------------------------------
for tool, parser in PARSERS.items():
    print(f"Parsing {tool} ...")
    ko_sets = {}   # sample → set of KOs
    for sample in samples:
        try:
            ko_sets[sample] = parser(sample)
        except Exception as e:
            warnings.warn(f"  [{tool}] {sample}: {e}")
            ko_sets[sample] = set()
        n = len(ko_sets[sample])
        print(f"  {sample}: {n} KOs")

    # Union of all KOs this tool ever detected, restricted to bacterial module KOs
    raw_kos = set().union(*ko_sets.values())
    all_kos = sorted(raw_kos & BACTERIAL_MODULE_KOS)
    n_dropped = len(raw_kos) - len(all_kos)
    print(f"  Total unique KOs across all samples: {len(raw_kos):,} raw  →  "
          f"{len(all_kos):,} in bacterial modules  ({n_dropped:,} outside module universe dropped)")

    # Build binary matrix: rows = samples, columns = KOs
    matrix = pd.DataFrame(
        index=pd.Index(samples, name="sample"),
        columns=pd.Index(all_kos, name="ko_id"),
        dtype=np.int8,
    )
    for sample, ko_set in ko_sets.items():
        matrix.loc[sample] = 0
        if ko_set:
            present = [k for k in ko_set if k in matrix.columns]
            matrix.loc[sample, present] = 1

    out_path = OUT_DIR / f"{tool}_binary.csv"
    matrix.to_csv(out_path)
    print(f"  Saved: {out_path}  (shape: {matrix.shape})")
    print()

print("All done. KO_MATRICES/ is ready for plot_ko_figures.py")
