#!/usr/bin/env python3
"""
BLIMMP wrapper for puma — uses puma-specific reference data paths.

Usage (same CLI as BLIMMP_Version_26Sep2026.py):
    python blimmp_puma_wrapper.py <domtblout> -f domtblout -s 0.99 \
        -o <output_prefix> [--verbose] [--no-substitutes]
"""

import sys
import os
import argparse
import logging
import importlib.util
from pathlib import Path

# ---------------------------------------------------------------------------
# BLIMMP source and reference data location on puma
# ---------------------------------------------------------------------------
BLIMMP_DIR = Path("/xdisk/twheeler/nsontakke/Removal_Study_BLIMMP/BLIMMP_8Sep2026/BLIMMP")
BLIMMP_SCRIPT = BLIMMP_DIR / "BLIMMP_Version_26Sep2026.py"

if not BLIMMP_SCRIPT.is_file():
    sys.exit(f"ERROR: BLIMMP script not found at {BLIMMP_SCRIPT}")

# Import BLIMMP as a module (the if __name__=='__main__' guard prevents
# main() from running automatically on import)
spec = importlib.util.spec_from_file_location("blimmp", BLIMMP_SCRIPT)
blimmp = importlib.util.module_from_spec(spec)
spec.loader.exec_module(blimmp)

# ---------------------------------------------------------------------------
# Reference data paths — all flat inside BLIMMP_DIR on puma
# ---------------------------------------------------------------------------
PATHS = blimmp.Paths(
    counts_dir               = BLIMMP_DIR / "ATB Frequency",
    onehop_dir               = BLIMMP_DIR / "ONE_HOP_NEIGHBOR_DATA-2",
    twohop_dir               = BLIMMP_DIR / "TWO_HOP_NEIGHBOR_DATA-2",
    module_neighbor_dir      = BLIMMP_DIR / "MODULE_ALL_NEIGHBOR_DATA",
    module_eq_json           = BLIMMP_DIR / "KEGG_Module_Equations_18AUG26.json",
    module_json_dir          = BLIMMP_DIR / "KEGG_Graphs_Generated_13AUG2026",
    kofam_ko_list_path       = BLIMMP_DIR / "ko_list.txt",
    module_frequencies       = BLIMMP_DIR / "module_freq.txt",
    module_reaction_dir      = BLIMMP_DIR / "module_ko_reaction.json",
    module_descriptions_path = BLIMMP_DIR / "kegg_bacteria_modules.json",
    ko_reaction_path         = BLIMMP_DIR / "ko_reaction.list",
    cooccurrence_lookup_dir  = BLIMMP_DIR,  # CoOccurrence_AllVsAll_domain_level_priors.json is here
)


def main():
    p = argparse.ArgumentParser(
        description="BLIMMP (puma) — metabolic pathway completeness via Bayesian graph inference.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("file",        help="Path to the HMMER .domtblout file")
    p.add_argument("-f", "--format", choices=["domtblout"], required=True)
    p.add_argument("-s", "--sigma",  type=float, required=True,
                   help="Genome completeness (0.0-1.0); use 0.99 for unknown/complete")
    p.add_argument("-t", "--taxonomy", default="bacteria", metavar="NAME",
                   help="Taxonomic group for priors (default: bacteria)")
    p.add_argument("-o", "--output", required=True, metavar="PREFIX",
                   help="Output prefix, e.g. BLIMMP_SPLICED/MED4_50pct_rep1/MED4_50pct_rep1_")
    p.add_argument("-v", "--verbose", action="store_true",
                   help="Write debug log + intermediate CSVs; outputs BOTH pass-1 and pass-2")
    p.add_argument("--no-substitutes", action="store_true",
                   help="Pass-1 only; skip substitution pipeline")
    args = p.parse_args()

    logfile_path = f"{args.output}_debug.log" if args.verbose else None
    logging.getLogger("numba").setLevel(logging.ERROR)

    logger = logging.getLogger()
    logger.handlers.clear()

    if args.verbose:
        logger.setLevel(logging.DEBUG)
        fh = logging.FileHandler(logfile_path, mode="w")
        fh.setLevel(logging.DEBUG)
        ch = logging.StreamHandler()
        ch.setLevel(logging.INFO)
        fmt = logging.Formatter("[%(levelname)s] %(message)s")
        fh.setFormatter(fmt); ch.setFormatter(fmt)
        logger.addHandler(fh); logger.addHandler(ch)
    else:
        logger.setLevel(logging.INFO)
        ch = logging.StreamHandler()
        ch.setLevel(logging.INFO)
        ch.setFormatter(logging.Formatter("[%(levelname)s] %(message)s"))
        logger.addHandler(ch)

    cfg = blimmp.RunConfig(
        input_file     = args.file,
        fmt            = args.format,
        sigma          = args.sigma,
        taxonomy       = args.taxonomy,
        output_prefix  = args.output,
        verbose        = args.verbose,
        logfile_path   = logfile_path,
        no_substitutes = args.no_substitutes,
    )

    blimmp.BlimmpPipeline(cfg, PATHS).run()


if __name__ == "__main__":
    main()