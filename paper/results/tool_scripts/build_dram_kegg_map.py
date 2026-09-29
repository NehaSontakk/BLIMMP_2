#!/usr/bin/env python3
"""
build_dram_kegg_map.py

Builds a DRAM module name → KEGG module ID (M#####) map directly from
DRAM's own module_step_form TSV — no text mining or KEGG API needed.

module_step_form columns used:
  module       → KEGG module ID  (M#####)
  module_name  → the exact English name DRAM uses in metabolism_summary.xlsx

Any DRAM module description that doesn't appear in module_step_form
(CRISPR subtypes, ribosomes, polymerases — KEGG BRITE entries) has no
M##### equivalent and is marked unmatched.

Usage:
  python3 build_dram_kegg_map.py \
      --form  /path/to/DRAM_data/module_step_form.20260727.tsv \
      --dram  /path/to/distill/metabolism_summary.xlsx \
      --out   dram_to_kegg_map.csv
"""

import argparse
import sys
import pandas as pd


def load_form(form_path: str) -> dict[str, str]:
    """Return {module_name → kegg_id} from module_step_form.tsv."""
    print(f"Reading form: {form_path}")
    df = pd.read_csv(form_path, sep="\t", usecols=["module", "module_name"],
                     dtype=str)
    df = df.dropna(subset=["module", "module_name"])
    df = df[df["module"].str.match(r"^M\d{5}$")]
    mapping = (df.drop_duplicates(subset=["module_name"])
                 .set_index("module_name")["module"]
                 .to_dict())
    print(f"  {len(mapping)} unique (module_name → M#####) pairs")
    return mapping


def load_dram_modules(xlsx_path: str) -> list[str]:
    print(f"Reading DRAM output: {xlsx_path}")
    df = pd.read_excel(xlsx_path, sheet_name="MISC")
    if "module" not in df.columns:
        sys.exit("ERROR: 'module' column not found in MISC sheet.")
    mods = df["module"].dropna().unique().tolist()
    print(f"  {len(mods)} unique module descriptions")
    return mods


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--form", required=True,
                    help="Path to DRAM_data/module_step_form.*.tsv")
    ap.add_argument("--dram", required=True,
                    help="Path to DRAM distill/metabolism_summary.xlsx")
    ap.add_argument("--out", default="dram_to_kegg_map.csv",
                    help="Output CSV (default: dram_to_kegg_map.csv)")
    args = ap.parse_args()

    form_map  = load_form(args.form)
    dram_mods = load_dram_modules(args.dram)

    rows = []
    n_matched = n_unmatched = 0
    for raw in dram_mods:
        kegg_id = form_map.get(raw, "")
        matched = bool(kegg_id)
        rows.append({"dram_module": raw, "kegg_id": kegg_id,
                     "matched": matched})
        if matched:
            n_matched += 1
        else:
            n_unmatched += 1

    print(f"\nMatched (have M#####):   {n_matched}")
    print(f"Unmatched (BRITE/other): {n_unmatched}")

    df_out = pd.DataFrame(rows)
    df_out.to_csv(args.out, index=False)
    print(f"Saved → {args.out}")

    unmatched = df_out[~df_out["matched"]]["dram_module"].tolist()
    if unmatched:
        print("\nNo M##### equivalent (KEGG BRITE / DRAM-specific):")
        for desc in unmatched:
            print(f"  {desc!r}")


if __name__ == "__main__":
    main()
