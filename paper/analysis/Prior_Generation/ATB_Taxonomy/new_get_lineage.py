# get_lineage.py
import os
import glob
import subprocess
import pandas as pd
import re
from typing import List, Optional

# ---------- CONFIG ----------
TAXONKIT_DATA = "/xdisk/twheeler/nsontakke/ATB_Taxonomy"
TANTAN_DIR = "/xdisk/twheeler/nsontakke/ATB_Analysis_0725/ATB_TANTAN"
MASTER_LIST = "/xdisk/twheeler/nsontakke/ATB_Analysis_0725/DOWNLOAD_SCRIPTS/file_list.all.20240805.tsv"
OUTPUT_A = "ATB_sample_species.tsv"                 # cleaned names
OUTPUT_B = "ATB_sample_species1.tsv"                # with taxid
OUTPUT_C = "ATB_sample_species_full_lineage.tsv"
CLEANED_OUTPUT = "ATB_sample_species_full_lineage_clean.tsv"
CLEAN_MAP = "ATB_species_clean_map.tsv"

# Strip genus/species “partitions” like _A, _B, _G3, _K, etc.
# Also strip on species epithets like "botulinum_A" or "thuringiensis_S".
_PARTED_RE = re.compile(r"^([A-Z][a-z]+)_[A-Za-z0-9]+$")      # for Genus
_PARTED_SPECIES_RE = re.compile(r"^([a-z][a-z\-]+)_[A-Za-z0-9]+$")  # for species epithet
_SP_NUM_RE = re.compile(r"^sp\d{3,}$")                        # sp followed by digits (3+)

def _strip_quotes(s: str) -> str:
    return str(s).strip().strip('"').strip("'")

def clean_taxon_name(name: str, prefer_genus_sp: bool = True) -> str:
    """
    - Remove quoted wrappers
    - Normalize 'Candidatus' / 'Ca.' prefix handling
    - Strip parted genus/species suffixes (_A, _B, _K, _S, _G3, etc.)
    - Replace tokens like sp010998615 -> 'sp.' (or drop to genus if prefer_genus_sp=False)
    """
    if name is None or not str(name).strip():
        return name
    s = _strip_quotes(name)
    parts = s.split()
    if not parts:
        return s

    # Handle Candidatus prefix (kept as-is, but genus is next token)
    cand_prefix = None
    genus_idx = 0
    if parts[0] in ("Candidatus", "Ca."):
        cand_prefix = parts[0]
        genus_idx = 1

    # Normalize genus (strip partitions)
    if len(parts) > genus_idx:
        genus = parts[genus_idx]
        m = _PARTED_RE.match(genus)
        if m:
            parts[genus_idx] = m.group(1)

    # Normalize species epithet if present
    species_idx = genus_idx + 1
    if len(parts) > species_idx:
        sp = parts[species_idx]

        # Replace spNNN… placeholders
        if _SP_NUM_RE.match(sp):
            if prefer_genus_sp:
                parts = parts[:species_idx] + ["sp."]  # "Genus sp."
            else:
                parts = parts[:species_idx]            # just "Genus"
        else:
            # Strip partition suffix on species epithet (e.g., thuringiensis_S)
            m2 = _PARTED_SPECIES_RE.match(sp)
            if m2:
                parts[species_idx] = m2.group(1)

    # If there is a third token (subspecies/strain like "strainX"), leave it alone.
    return " ".join(parts)

def run_cmd(cmd: List[str], input_text: Optional[str] = None) -> str:
    res = subprocess.run(
        cmd,
        input=input_text.encode() if input_text is not None else None,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if res.returncode != 0:
        raise RuntimeError(f"Command failed: {' '.join(cmd)}\nSTDERR:\n{res.stderr.decode()}")
    return res.stdout.decode()

def name2taxids(qname: str) -> List[str]:
    """
    Returns a list of candidate taxids for the given name using taxonkit.
    """
    if not qname or not qname.strip():
        return []
    cmd = [
        "taxonkit", "name2taxid",
        "--data-dir", TAXONKIT_DATA,
        "--name-field", "1"
    ]
    out = run_cmd(cmd, input_text=qname + "\n").strip()
    # Output format: "<name>\t<taxid[;taxid2;...]>"
    # Handle empty returns gracefully
    if not out:
        return []
    fields = out.split("\t")
    if len(fields) < 2 or not fields[1].strip():
        return []
    # split by ';' or ',' just in case
    ids = re.split(r"[;,]", fields[1].strip())
    return [tid.strip() for tid in ids if tid.strip().isdigit()]

def lca_of_taxids(taxids: List[str]) -> Optional[str]:
    """
    Compute LCA taxid using taxonkit lca.
    """
    if not taxids:
        return None
    if len(taxids) == 1:
        return taxids[0]
    # taxonkit lca expects one line with space-separated taxids
    line = " ".join(taxids) + "\n"
    cmd = ["taxonkit", "lca", "--data-dir", TAXONKIT_DATA]
    out = run_cmd(cmd, input_text=line).strip()
    # Output is "taxid\trank\tname" (single line), or may include the input as first column in some versions.
    if not out:
        return None
    cols = out.split("\t")
    # Heuristics: last numeric column is usually the LCA taxid
    for c in cols:
        if c.isdigit():
            return c
    return None

def main():
    # 1. Gather masked filenames and extract sample IDs
    mask_files = glob.glob(os.path.join(TANTAN_DIR, "*_masked.fa"))
    sample_ids = [os.path.basename(fp).replace("_masked.fa", "") for fp in mask_files]
    print(f"Found {len(sample_ids)} TANTAN samples")

    # 2. Load the master file list and filter by sample IDs
    df_all = pd.read_csv(MASTER_LIST, sep="\t", dtype=str)
    print(f"Loaded master file list: {len(df_all)} rows")
    df_all.columns = df_all.columns.str.strip()
    if 'sample' not in df_all.columns or 'species_sylph' not in df_all.columns:
        print("ERROR: Expected columns 'sample' and 'species_sylph' not found.")
        print("Available columns:", df_all.columns.tolist())
        return

    filtered = df_all[df_all['sample'].isin(sample_ids)].copy()
    print(f"Filtered to {len(filtered)} matching rows")

    # 3. Clean names and write audit map
    filtered['species_clean'] = filtered['species_sylph'].apply(
        lambda x: clean_taxon_name(x, prefer_genus_sp=True)
    )

    filtered[['species_sylph', 'species_clean']].drop_duplicates().to_csv(
        CLEAN_MAP, sep="\t", index=False
    )

    # 4. Resolve taxids with LCA fallback for ambiguous names
    taxids = []
    for nm in filtered['species_clean']:
        cands = name2taxids(nm)
        if not cands:
            # Fallback: try collapsing "Genus sp." to just Genus
            nm2 = re.sub(r"\s+sp\.\s*$", "", nm)
            cands = name2taxids(nm2) if nm2 != nm else []
        lca = lca_of_taxids(cands) if cands else ""
        taxids.append(lca if lca is not None else "")

    filtered['taxid'] = taxids

    # 5. Write cleaned sample/species for record
    filtered[['sample', 'species_clean']].to_csv(OUTPUT_A, sep="\t", index=False, header=False)
    print(f"Wrote sample/species pairs (cleaned) to {OUTPUT_A}")

    # 6. Write taxid-augmented file
    filtered[['sample', 'species_clean', 'taxid']].to_csv(OUTPUT_B, sep="\t", index=False, header=False)
    print(f"Wrote taxid-augmented file to {OUTPUT_B}")

    # 7. Lineage using taxid field
    cmd_lineage = (
        f"taxonkit lineage "
        f"--data-dir {TAXONKIT_DATA} "
        f"--taxid-field 3 "
        f"--show-lineage-ranks "
        f"--show-name "
        f"{OUTPUT_B} "
        f"> {OUTPUT_C}"
    )
    subprocess.run(cmd_lineage, shell=True, executable="/bin/bash", check=True)
    print(f"Wrote full lineage file to {OUTPUT_C}")

    # 8. Remove repeated header lines intermixed in the output
    with open(OUTPUT_C, 'r') as infile, open(CLEANED_OUTPUT, 'w') as outfile:
        for line in infile:
            if line.startswith("cellular root"):
                continue
            outfile.write(line)
    print(f"Removed repeated headers; cleaned file: {CLEANED_OUTPUT}")

    # 9. Preview first 10 lines
    print("\nFirst 10 lines of cleaned full lineage output:")
    with open(CLEANED_OUTPUT, 'r') as f:
        for i in range(10):
            ln = f.readline()
            if not ln:
                break
            print(ln.rstrip("\n"))

if __name__ == "__main__":
    main()
