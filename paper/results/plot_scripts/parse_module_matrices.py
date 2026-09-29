#!/usr/bin/env python3
"""
parse_module_matrices.py

Convert all tool outputs into binary MODULE presence/absence matrices.
Run ONCE (reads many files); outputs go to MODULE_MATRICES/.

Module presence rules per tool:
  anvi'o          — pathwise_module_is_complete == True (per sample, any path)
  BLIMMP          — module_confidence >= 0.60
                    (reads *__BLIMMP_substituted_module_probabilities.csv)
  MetaPathPredict — binary prediction >= 0.5
  METABOLIC       — KO→module stepwise completeness >= 0.75
  DRAM            — KO→module stepwise completeness >= 0.75
  KEMET           — KO→module stepwise completeness >= 0.75

All modules are restricted to the bacterial KEGG module universe in
module_ko_reaction.json (same universe BLIMMP searches).

Usage:
    conda activate blimmp-work
    python parse_module_matrices.py
"""

import json
import re
import warnings
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Paths & config
# ---------------------------------------------------------------------------
BASE_DIR = Path("/xdisk/twheeler/nsontakke/Removal_Study_BLIMMP")
OUT_DIR  = BASE_DIR / "MODULE_MATRICES"
OUT_DIR.mkdir(exist_ok=True)

MODULE_KO_JSON = BASE_DIR / "BLIMMP_8Sep2026" / "BLIMMP" / "module_ko_reaction.json"
with open(MODULE_KO_JSON) as f:
    MODULE_KO_DATA = json.load(f)

ALL_MODULES = sorted(MODULE_KO_DATA.keys())
MODULE_KOS  = {m: set(d.keys()) for m, d in MODULE_KO_DATA.items()}
print(f"Loaded {len(ALL_MODULES)} modules from module_ko_reaction.json\n")

BLIMMP_THRESH = 0.60   # module_confidence threshold for BLIMMP
KO_THRESH     = 0.75   # stepwise completeness threshold for KO-based tools

SAMPLE_RE = re.compile(r"^(.+)_(\d+)percentremoved_replicate(\d+)$")

def parse_sample_name(name):
    m = SAMPLE_RE.match(name)
    return (m.group(1), int(m.group(2)), int(m.group(3))) if m else None

def ko_set_to_modules(ko_set: set, threshold: float = KO_THRESH) -> dict:
    """Convert a detected KO set to binary module presence via stepwise completeness."""
    result = {}
    for mod_id in ALL_MODULES:
        mod_kos = MODULE_KOS[mod_id]
        if not mod_kos:
            result[mod_id] = 0
        else:
            comp = len(mod_kos & ko_set) / len(mod_kos)
            result[mod_id] = 1 if comp >= threshold else 0
    return result

# ---------------------------------------------------------------------------
# Discover samples from KOFAM_SPLICED (most complete sample listing)
# ---------------------------------------------------------------------------
print("Discovering samples from KOFAM_SPLICED ...")
kofam_root = BASE_DIR / "KOFAM_SPLICED"
all_dirs = sorted(d.name for d in kofam_root.iterdir() if d.is_dir())

samples, meta_rows = [], []
for name in all_dirs:
    parsed = parse_sample_name(name)
    if parsed is None:
        continue
    organism, pct, rep = parsed
    samples.append(name)
    meta_rows.append(dict(sample=name, organism=organism,
                          removal_pct=pct, replicate=rep))

meta_df = pd.DataFrame(meta_rows).set_index("sample")
meta_df.to_csv(OUT_DIR / "sample_metadata.csv")
print(f"  {len(samples)} samples  →  {OUT_DIR}/sample_metadata.csv\n")

# ---------------------------------------------------------------------------
# Per-tool parsers
# ---------------------------------------------------------------------------

def parse_anvio(sample: str) -> dict:
    """pathwise_module_is_complete (or stepwise fallback)."""
    p = BASE_DIR / "ANVIO_SPLICED" / sample / f"{sample}_modules.txt"
    result = {m: 0 for m in ALL_MODULES}
    if not p.exists():
        return result
    df = pd.read_csv(p, sep="\t", low_memory=False)
    if "pathwise_module_is_complete" in df.columns:
        col = "pathwise_module_is_complete"
    elif "stepwise_module_is_complete" in df.columns:
        col = "stepwise_module_is_complete"
        print(f"    [WARN] Using stepwise for {sample}")
    else:
        return result
    df["present"] = df[col].astype(str).str.strip() == "True"
    # Multiple paths per module → any path complete counts
    pres = df.groupby("module")["present"].any()
    for mod, v in pres.items():
        if str(mod) in result:
            result[str(mod)] = 1 if v else 0
    return result


def _load_module_probs(sample: str) -> Optional[pd.DataFrame]:
    """Load the PRE-SUBSTITUTION module probabilities file.
    Contains module_confidence (post-network) and module_confidence_prior (pre-network).

    IMPORTANT: explicitly excludes *substituted* files — the glob
    '*_module_probabilities.csv' would otherwise also match
    '*_substituted_module_probabilities.csv'.
    """
    d = BASE_DIR / "BLIMMP_SPLICED" / sample
    if not d.exists():
        return None
    for pat, sep in [("*_module_probabilities.csv", ","),
                     ("*_module_probabilities_v2.tsv", "\t"),
                     ("*_module_probabilities.tsv", "\t")]:
        # Exclude any file that has "substituted" in its name
        hits = [h for h in d.glob(pat) if "substituted" not in h.name]
        if hits:
            return pd.read_csv(hits[0], sep=sep)
    return None


def _apply_thresh(df: pd.DataFrame, col: str) -> dict:
    """Vectorised helper: returns {module: 0/1} for df[col] >= BLIMMP_THRESH."""
    result = {m: 0 for m in ALL_MODULES}
    sub = df[["module", col]].copy()
    sub["module"] = sub["module"].astype(str).str.strip()
    sub = sub[sub["module"].isin(result)]
    for mod in sub.loc[sub[col] >= BLIMMP_THRESH, "module"]:
        result[mod] = 1
    return result


def parse_blimmp_prior(sample: str) -> dict:
    """*_module_probabilities.csv → module_confidence_prior >= 0.60.
    Direct HMM evidence only — before neighbor network update."""
    df = _load_module_probs(sample)
    if df is None or "module_confidence_prior" not in df.columns:
        return {m: 0 for m in ALL_MODULES}
    return _apply_thresh(df, "module_confidence_prior")


def parse_blimmp_nosub(sample: str) -> dict:
    """*_module_probabilities.csv → module_confidence >= 0.60.
    After neighbor network update, before substitution."""
    df = _load_module_probs(sample)
    if df is None or "module_confidence" not in df.columns:
        return {m: 0 for m in ALL_MODULES}
    return _apply_thresh(df, "module_confidence")


def parse_blimmp(sample: str) -> dict:
    """*__BLIMMP_substituted_module_probabilities.csv → module_confidence >= 0.60.
    Full pipeline: direct HMM + neighbor update + substitution."""
    d = BASE_DIR / "BLIMMP_SPLICED" / sample
    result = {m: 0 for m in ALL_MODULES}
    if not d.exists():
        return result
    hits = list(d.glob("*__BLIMMP_substituted_module_probabilities.csv"))
    if not hits:
        return result
    df = pd.read_csv(hits[0], usecols=["module", "module_confidence"])
    return _apply_thresh(df, "module_confidence")


_MPP_CACHE = None

def parse_mpp(sample: str) -> dict:
    """Binary prediction >= 0.5 from metapathpredict_all.tsv."""
    global _MPP_CACHE
    result = {m: 0 for m in ALL_MODULES}
    if _MPP_CACHE is None:
        path = BASE_DIR / "METAPATHPREDICT_SPLICED" / "metapathpredict_all.tsv"
        if not path.exists():
            print(f"  [MISSING] MetaPathPredict: {path}")
            _MPP_CACHE = pd.DataFrame()
            return result
        df = pd.read_csv(path, sep="\t", index_col=0)
        new_idx = []
        for idx in df.index:
            parts = str(idx).replace("\\", "/").split("/")
            new_idx.append(parts[-2] if len(parts) >= 2 else parts[-1])
        df.index = new_idx
        print(f"  MetaPathPredict: {df.shape[0]} samples × {df.shape[1]} modules")
        _MPP_CACHE = df
    if _MPP_CACHE.empty or sample not in _MPP_CACHE.index:
        return result
    row = _MPP_CACHE.loc[sample].astype(float)
    for mod_id in row.index:
        if str(mod_id) in result:
            result[str(mod_id)] = 1 if row[mod_id] >= 0.5 else 0
    return result


def parse_metabolic_kos(sample: str) -> set:
    """total.hits.txt: col[0] = KO when line has >= 2 space-separated tokens."""
    p = (BASE_DIR / "METABOLIC_SPLICED"
         / f"METABOLIC_{sample}" / "KEGG_identifier_result" / "total.hits.txt")
    kos = set()
    if not p.exists():
        return kos
    with open(p) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                kos.add(parts[0].strip())
    return kos


def parse_dram_kos(sample: str) -> set:
    """annotations.tsv: ko_id column."""
    p = BASE_DIR / "DRAM_SPLICED" / sample / "annotations.tsv"
    if not p.exists():
        return set()
    df = pd.read_csv(p, sep="\t", usecols=["ko_id"])
    return {k.strip() for k in df["ko_id"].dropna().astype(str)
            if re.match(r"^K\d{5}$", k.strip())}


def parse_kemet_kos(sample: str) -> set:
    """ktests/*.ktest: one KO per line."""
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


# ---------------------------------------------------------------------------
# Tool configs
# ---------------------------------------------------------------------------
# Tools whose output is already module-level
DIRECT_TOOLS = {
    "anvio":           parse_anvio,
    "BLIMMP_prior":    parse_blimmp_prior,    # direct HMM evidence only
    "BLIMMP_nosub":    parse_blimmp_nosub,    # + neighbor network update
    "BLIMMP":          parse_blimmp,          # + substitution (full pipeline)
    "MetaPathPredict": parse_mpp,
}

# Tools that output KO sets → converted to modules via stepwise completeness >= KO_THRESH
KO_TOOLS = {
    "METABOLIC": parse_metabolic_kos,
    "DRAM":      parse_dram_kos,
    "KEMET":     parse_kemet_kos,
}

# ---------------------------------------------------------------------------
# Build and save matrices
# ---------------------------------------------------------------------------
def build_and_save(tool: str, sample_results: dict):
    df = pd.DataFrame.from_dict(
        sample_results, orient="index", columns=ALL_MODULES
    ).fillna(0).astype(np.int8)
    df.index.name = "sample"
    out_path = OUT_DIR / f"{tool}_module_binary.csv"
    df.to_csv(out_path)
    n_nonzero = (df.sum(axis=1) > 0).sum()
    print(f"  Saved: {out_path}  shape={df.shape}"
          f"  ({n_nonzero}/{len(df)} samples have ≥1 module)")


for tool, parser in DIRECT_TOOLS.items():
    print(f"Parsing {tool} (direct module output) ...")
    sample_results = {}
    for sample in samples:
        try:
            sample_results[sample] = parser(sample)
        except Exception as e:
            warnings.warn(f"  [{tool}] {sample}: {e}")
            sample_results[sample] = {m: 0 for m in ALL_MODULES}
    build_and_save(tool, sample_results)
    print()

for tool, ko_parser in KO_TOOLS.items():
    print(f"Parsing {tool} (KO→module, stepwise completeness ≥{KO_THRESH}) ...")
    sample_results = {}
    for sample in samples:
        try:
            ko_set = ko_parser(sample)
            sample_results[sample] = ko_set_to_modules(ko_set, KO_THRESH)
        except Exception as e:
            warnings.warn(f"  [{tool}] {sample}: {e}")
            sample_results[sample] = {m: 0 for m in ALL_MODULES}
    build_and_save(tool, sample_results)
    print()

print("Done. MODULE_MATRICES/ is ready for plot_module_recovery.py / plot_module_matrix.py")