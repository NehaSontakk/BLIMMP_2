#!/usr/bin/env python3
import csv, json, glob
from pathlib import Path

BASE = Path("/xdisk/cgoubert/nsontakke/ATB_KO_Frequencies")

INTERSECT_DIR = BASE / "Lineage_Specific_Data" / "KO_Intersections_Lineage_Specific"
COUNTS_DIR    = BASE / "Lineage_Specific_Data" / "KO_Counts_Lineage_Specific"
NEIGHBORS_JSON = BASE / "One_Hop_Neighbor_Adjacency.json"

OUT_PI_GIVEN_J = BASE / "Lineage_Specific_Data" / "KO_Conditionals_Pi_given_j"
OUT_PJ_GIVEN_I = BASE / "Lineage_Specific_Data" / "KO_Conditionals_Pj_given_i"
OUT_PI_GIVEN_J.mkdir(parents=True, exist_ok=True)
OUT_PJ_GIVEN_I.mkdir(parents=True, exist_ok=True)


def base_from_intersection(p: Path) -> str:
    # ko_intersection_ko_matrix_<basename>.tsv -> ko_matrix_<basename>
    return p.stem[len("ko_intersection_"):]


def load_counts(path: Path) -> dict:
    """Read KO counts file (2 cols: KOs, KO_count) into dict."""
    counts = {}
    with path.open() as f:
        rdr = csv.reader(f, delimiter="\t")
        header = next(rdr)
        # find indexes
        try:
            idx_ko = header.index("KOs")
        except ValueError:
            idx_ko = 0
        try:
            idx_cnt = header.index("KO_count")
        except ValueError:
            idx_cnt = 1
        for row in rdr:
            if not row:
                continue
            counts[row[idx_ko]] = int(row[idx_cnt])
    return counts


def safe_div(n, d):
    return (float(n) / d) if d else 0.0


# ---- load neighbor dictionary ----
with open(NEIGHBORS_JSON) as f:
    NEIGHBORS = json.load(f)

# ---- iterate intersection files ----
for f_inter in sorted(glob.glob(str(INTERSECT_DIR / "ko_intersection_ko_matrix_*.tsv"))):
    f_inter = Path(f_inter)
    base = base_from_intersection(f_inter)

    counts_path = COUNTS_DIR / f"ko_counts_{base}.tsv"
    if not counts_path.exists():
        print(f"[WARN] missing counts for {base}")
        continue
    counts = load_counts(counts_path)

    # open outputs
    out_pi = (OUT_PI_GIVEN_J / f"neighbors_Pi_given_j_{base}.tsv").open("w", newline="")
    out_pj = (OUT_PJ_GIVEN_I / f"neighbors_Pj_given_i_{base}.tsv").open("w", newline="")
    w_pi = csv.writer(out_pi, delimiter="\t")
    w_pj = csv.writer(out_pj, delimiter="\t")
    w_pi.writerow(["KO_j","KO_i","co_ij","count_j","P_i_given_j"])
    w_pj.writerow(["KO_i","KO_j","co_ij","count_i","P_j_given_i"])

    with f_inter.open(newline="") as fin:
        rdr = csv.reader(fin, delimiter="\t")
        header = next(rdr)
        col_to_idx = {k:i for i,k in enumerate(header)}
        present = set(header)

        # only consider neighbors that are present
        valid_j = [k for k in NEIGHBORS if k in present]
        # map all needed neighbor columns
        needed_cols = {}
        for j in valid_j:
            for i in NEIGHBORS[j]:
                if i in present:
                    needed_cols[i] = col_to_idx[i]

        for row in rdr:
            if not row: continue
            KO_row = row[0]
            if KO_row not in valid_j:
                continue
            cnt_j = counts.get(KO_row, 0)
            for KO_i in NEIGHBORS.get(KO_row, []):
                col_idx = needed_cols.get(KO_i)
                if col_idx is None: 
                    continue
                try:
                    co_ij = int(row[col_idx]) if row[col_idx] else 0
                except ValueError:
                    co_ij = 0
                cnt_i = counts.get(KO_i, 0)

                w_pi.writerow([KO_row, KO_i, co_ij, cnt_j, f"{safe_div(co_ij,cnt_j):.10f}"])
                w_pj.writerow([KO_i, KO_row, co_ij, cnt_i, f"{safe_div(co_ij,cnt_i):.10f}"])

    out_pi.close()
    out_pj.close()
    print(f"[OK] {f_inter.name} -> wrote {out_pi.name}, {out_pj.name}")
