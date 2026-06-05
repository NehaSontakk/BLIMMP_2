# new_get_lineage_fast.py
import os
import glob
import subprocess
import pandas as pd
import re
from typing import List, Dict, Optional

# ---------- CONFIG ----------
TAXONKIT_DATA = "/xdisk/twheeler/nsontakke/ATB_Taxonomy"
TANTAN_DIR = "/xdisk/twheeler/nsontakke/ATB_Analysis_0725/ATB_TANTAN"
MASTER_LIST = "/xdisk/twheeler/nsontakke/ATB_Analysis_0725/DOWNLOAD_SCRIPTS/file_list.all.20240805.tsv"
OUTPUT_A = "ATB_sample_species.tsv"                  # cleaned names
OUTPUT_B = "ATB_sample_species1.tsv"                 # with taxid
OUTPUT_C = "ATB_sample_species_full_lineage.tsv"
CLEANED_OUTPUT = "ATB_sample_species_full_lineage_clean.tsv"
CLEAN_MAP = "ATB_species_clean_map.tsv"

# threads for taxonkit
THREADS = str(max(1, (os.cpu_count() or 4) // 2))

# ---------- Cleaning rules ----------
_PARTED_GENUS_RE = re.compile(r"^([A-Z][a-z]+)_[A-Za-z0-9]+$")
_PARTED_SPECIES_RE = re.compile(r"^([a-z][a-z\-]+)_[A-Za-z0-9]+$")
_SP_NUM_RE = re.compile(r"^sp\d{3,}$")

def _strip_quotes(s: str) -> str:
    return str(s).strip().strip('"').strip("'")

def clean_taxon_name(name: str, prefer_genus_sp: bool = True) -> str:
    """
    Strip genus/species partitions (_A, _B, _K, _S, _G3...), fold sp######## -> 'sp.' (or drop to genus).
    Handles 'Candidatus' / 'Ca.' prefix.
    """
    if name is None or not str(name).strip():
        return ""
    s = _strip_quotes(name)
    parts = s.split()
    if not parts:
        return s

    cand_prefix = None
    genus_idx = 0
    if parts[0] in ("Candidatus", "Ca."):
        cand_prefix = parts[0]
        genus_idx = 1

    # genus
    if len(parts) > genus_idx:
        g = parts[genus_idx]
        m = _PARTED_GENUS_RE.match(g)
        if m:
            parts[genus_idx] = m.group(1)

    # species
    sp_idx = genus_idx + 1
    if len(parts) > sp_idx:
        sp = parts[sp_idx]
        if _SP_NUM_RE.match(sp):
            if prefer_genus_sp:
                parts = parts[:sp_idx] + ["sp."]
            else:
                parts = parts[:sp_idx]
        else:
            m2 = _PARTED_SPECIES_RE.match(sp)
            if m2:
                parts[sp_idx] = m2.group(1)

    return " ".join(parts)

# ---------- Helpers ----------
def run_cmd(cmd: List[str], input_text: Optional[str] = None) -> str:
    res = subprocess.run(
        cmd,
        input=(input_text.encode() if input_text is not None else None),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if res.returncode != 0:
        raise RuntimeError(
            f"Command failed: {' '.join(cmd)}\nSTDERR:\n{res.stderr.decode()}"
        )
    return res.stdout.decode()

def name2taxids_bulk(names: List[str]) -> Dict[str, List[str]]:
    """
    Bulk resolve names to taxids using:
      taxonkit name2taxid --data-dir ... -i 1 -j THREADS
    Input: one name per line.
    Output: map name -> [taxid, taxid2, ...]  (empty list if not found)
    """
    if not names:
        return {}
    # TaxonKit reads tab-separated; but with -i 1, a single field line is fine.
    cmd = [
        "taxonkit", "name2taxid",
        "--data-dir", TAXONKIT_DATA,
        "--name-field", "1",
        "-j", THREADS,
    ]
    input_text = "\n".join(names) + "\n"
    out = run_cmd(cmd, input_text=input_text).strip().splitlines()

    mapping: Dict[str, List[str]] = {}
    for line in out:
        # Expect "<name>\t<taxids>" or "<name>\t"
        parts = line.split("\t")
        if not parts:
            continue
        nm = parts[0].strip()
        if len(parts) < 2 or not parts[1].strip():
            mapping[nm] = []
        else:
            ids = re.split(r"[;,]", parts[1].strip())
            mapping[nm] = [tid for tid in (i.strip() for i in ids) if tid.isdigit()]
    # Ensure all input names present
    for nm in names:
        mapping.setdefault(nm, [])
    return mapping

def lca_bulk(list_of_taxid_lists: List[List[str]]) -> List[str]:
    """
    Compute LCAs for many sets of taxids at once using:
      taxonkit lca -j THREADS
    Each input row: "tid1 tid2 ..." (or single tid).
    For empty lists: return "".
    """
    if not list_of_taxid_lists:
        return []
    lines = []
    for tids in list_of_taxid_lists:
        if not tids:
            lines.append("")  # placeholder
        else:
            lines.append(" ".join(tids))
    # Filter out empties before calling lca; keep index to restore order
    idx_map = [i for i, L in enumerate(list_of_taxid_lists) if L]
    nonempty_lines = [lines[i] for i in idx_map]
    lca_results = [""] * len(list_of_taxid_lists)

    if nonempty_lines:
        cmd = ["taxonkit", "lca", "--data-dir", TAXONKIT_DATA, "-j", THREADS]
        out = run_cmd(cmd, input_text="\n".join(nonempty_lines) + "\n").strip().splitlines()
        # taxonkit lca outputs: "taxid\trank\tname" (one per line)
        # We take the first numeric field on each line as the LCA taxid.
        for k, line in enumerate(out):
            lca_tid = ""
            for field in line.split("\t"):
                if field.isdigit():
                    lca_tid = field
                    break
            lca_results[idx_map[k]] = lca_tid

    return lca_results

def main():
    # 1) Gather masked filenames and sample IDs
    mask_files = glob.glob(os.path.join(TANTAN_DIR, "*_masked.fa"))
    sample_ids = [os.path.basename(fp).replace("_masked.fa", "") for fp in mask_files]
    print(f"Found {len(sample_ids)} TANTAN samples")

    # 2) Load only needed columns from the big TSV
    usecols = ["sample", "species_sylph"]
    df_all = pd.read_csv(MASTER_LIST, sep="\t", dtype=str, usecols=usecols)
    df_all.columns = df_all.columns.str.strip()
    print(f"Loaded master file list: {len(df_all)} rows")

    filtered = df_all[df_all['sample'].isin(sample_ids)].copy()
    print(f"Filtered to {len(filtered)} matching rows")

    # 3) Clean names
    filtered['species_clean'] = filtered['species_sylph'].apply(
        lambda x: clean_taxon_name(x, prefer_genus_sp=True)
    )
    filtered[['species_sylph', 'species_clean']].drop_duplicates().to_csv(
        CLEAN_MAP, sep="\t", index=False
    )

    # 4) Bulk name2taxid on unique cleaned names
    unique_names = sorted(set(filtered['species_clean']))
    print(f"Unique cleaned names: {len(unique_names)}")

    nm2taxids = name2taxids_bulk(unique_names)

    # 5) Bulk LCA for ambiguous names; pass singletons as well (cheap)
    taxid_lists = [nm2taxids[nm] for nm in unique_names]
    lcas = lca_bulk(taxid_lists)
    name_to_lca = {nm: lca for nm, lca in zip(unique_names, lcas)}

    # 6) Attach taxid to rows
    filtered['taxid'] = filtered['species_clean'].map(name_to_lca).fillna("")

    # 7) Write cleaned sample/species
    filtered[['sample', 'species_clean']].to_csv(OUTPUT_A, sep="\t", index=False, header=False)
    print(f"Wrote sample/species pairs (cleaned) to {OUTPUT_A}")

    # 8) Write taxid-augmented file
    filtered[['sample', 'species_clean', 'taxid']].to_csv(OUTPUT_B, sep="\t", index=False, header=False)
    print(f"Wrote taxid-augmented file to {OUTPUT_B}")

    # 9) Lineage using taxid field (bulk, multithreaded)
    cmd_lineage = (
        f"taxonkit lineage "
        f"--data-dir {TAXONKIT_DATA} "
        f"--taxid-field 3 "
        f"--show-lineage-ranks "
        f"-j {THREADS} "
        f"{OUTPUT_B} "
        f"> {OUTPUT_C}"
    )
    subprocess.run(cmd_lineage, shell=True, executable="/bin/bash", check=True)
    print(f"Wrote full lineage file to {OUTPUT_C}")

    # 10) Remove repeated header lines
    with open(OUTPUT_C, 'r') as infile, open(CLEANED_OUTPUT, 'w') as outfile:
        for line in infile:
            if line.startswith("cellular root"):
                continue
            outfile.write(line)
    print(f"Removed repeated headers; cleaned file: {CLEANED_OUTPUT}")

    # 11) Preview
    print("\nFirst 10 lines of cleaned full lineage output:")
    with open(CLEANED_OUTPUT, 'r') as f:
        for i in range(10):
            ln = f.readline()
            if not ln:
                break
            print(ln.rstrip("\n"))

if __name__ == "__main__":
    main()
