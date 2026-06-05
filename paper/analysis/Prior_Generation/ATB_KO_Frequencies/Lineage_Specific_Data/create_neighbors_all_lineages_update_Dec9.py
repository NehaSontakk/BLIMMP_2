#!/usr/bin/env python3
# Batch-refill One-Hop and Two-Hop neighbor JSONs for ALL lineage-specific files
# + Debug printouts for (row K00150, col K15634) and (row K11389, col K15634)
# + Output filenames strip "ko_intersection_ko_matrix_sampleids_" prefix

import json, csv, re
from pathlib import Path
from collections import defaultdict

BASE = Path("/xdisk/twheeler/nsontakke/copied_from_cgoubert/ATB_KO_Frequencies/Lineage_Specific_Data")

SAMPLE_ONE_HOP_JSON = BASE / "Sample_One_Hop_Neighbor.txt"
SAMPLE_TWO_HOP_JSON = BASE / "Sample_Two_Hop_Neighbor.txt"
SAMPLE_ALL_JSON = BASE / "Module_Neighborhood/Module_Based_Neighbor.txt"


COUNTS_DIR     = BASE / "KO_Counts_Lineage_Specific"
INTERSECT_DIR  = BASE / "KO_Intersections_Lineage_Specific"

OUT_ONE_DIR    = BASE / "ONE_HOP_NEIGHBOR_DATA"
OUT_TWO_DIR    = BASE / "TWO_HOP_NEIGHBOR_DATA"
OUT_ALL_DIR    = BASE / "MODULE_ALL_NEIGHBOR_DATA"

OUT_ONE_DIR.mkdir(parents=True, exist_ok=True)
OUT_TWO_DIR.mkdir(parents=True, exist_ok=True)
OUT_ALL_DIR.mkdir(parents=True, exist_ok=True)
# ---------------------------------------------------------------

ko_token = re.compile(r'^(K\d{5})(?:_\d+)?$') 

def base_ko(x: str) -> str | None:
    if x is None:
        return None
    m = ko_token.match(str(x).strip())
    return m.group(1) if m else None

def load_counts(path: Path) -> dict[str, float]:
    counts = {}
    if not path.exists():
        return counts
    with open(path, "r", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        for row in reader:
            if not row or len(row) < 2:
                continue
            ko = base_ko(row[0])
            if not ko:
                continue
            try:
                counts[ko] = float(row[1])
            except Exception:
                pass
    return counts

def collect_needed_pairs_from_samples(*blocks: dict) -> set[tuple[str,str]]:
    pairs = set()
    for block in blocks:
        for src, nb_dict in block.items():
            src_b = base_ko(src)
            if not src_b or not isinstance(nb_dict, dict):
                continue
            for nb in nb_dict.keys():
                if nb == "_count":
                    continue
                nb_b = base_ko(nb)
                if not nb_b or src_b == nb_b:
                    continue
                pairs.add(tuple(sorted((src_b, nb_b))))
    return pairs


def load_intersections_sparse(path_tsv: Path, needed_pairs: set[tuple[str,str]]) -> dict[tuple[str,str], float]:
    inter = {}
    if not path_tsv.exists() or not needed_pairs:
        return inter

    need_by_row = defaultdict(set)
    needed_cols = set()
    for a, b in needed_pairs:
        need_by_row[a].add(b)
        need_by_row[b].add(a)
        needed_cols.add(a); needed_cols.add(b)

    with open(path_tsv, "r", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader, None)
        if not header or len(header) < 2:
            return inter
        # normalize header KOs once
        raw_cols = header[1:]
        col_kos = [base_ko(x) for x in raw_cols]
        col_idx = {ko: (i+1) for i, ko in enumerate(col_kos) if ko and ko in needed_cols}

        for row in reader:
            if not row:
                continue
            row_ko = base_ko(row[0])
            if not row_ko:
                continue
            partners = need_by_row.get(row_ko)
            if not partners:
                continue
            for partner in partners:
                j = col_idx.get(partner)
                if j is None or j >= len(row):
                    continue
                cell = row[j]
                if not cell:
                    continue
                try:
                    v = float(cell)
                except ValueError:
                    continue
                inter[tuple(sorted((row_ko, partner)))] = v
    return inter

# --- direct, single-value fetch from intersection matrix for debugging ---
def get_matrix_value(path_tsv: Path, row_ko: str, col_ko: str):
    """Return float at (row=row_ko, col=col_ko) from wide TSV; None if missing."""
    if not path_tsv.exists():
        return None
    want_row = base_ko(row_ko)
    want_col = base_ko(col_ko)
    if not want_row or not want_col:
        return None

    with open(path_tsv, "r", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader, None)
        if not header or len(header) < 2:
            return None
        # map first matching normalized column
        col_kos = [base_ko(x) for x in header[1:]]
        col_idx = None
        for i, kk in enumerate(col_kos, start=1):
            if kk == want_col:
                col_idx = i
                break
        if col_idx is None:
            return None

        for row in reader:
            if not row:
                continue
            rk = base_ko(row[0])
            if rk != want_row:
                continue
            if col_idx >= len(row):
                return None
            cell = row[col_idx]
            try:
                return float(cell)
            except Exception:
                return None
    return None

def refill_counts(block: dict, counts: dict[str,float]) -> dict:
    out = {}
    for src, nb_dict in block.items():
        if not isinstance(nb_dict, dict):
            continue
        src_b = base_ko(src)
        refreshed = dict(nb_dict)
        refreshed["_count"] = float(counts.get(src_b, 0.0)) if src_b else float(refreshed.get("_count", 0.0))
        out[src] = refreshed
    return out

def refill_weights(block: dict, intersections: dict[tuple[str,str], float]) -> dict:
    out = {}
    for src, nb_dict in block.items():
        if not isinstance(nb_dict, dict):
            continue
        src_b = base_ko(src)
        if not src_b:
            out[src] = dict(nb_dict)
            continue
        refreshed = {}
        if "_count" in nb_dict:
            refreshed["_count"] = nb_dict["_count"]
        for nb, _old in nb_dict.items():
            if nb == "_count":
                continue
            nb_b = base_ko(nb)
            if not nb_b:
                refreshed[nb] = 0.0
                continue
            refreshed[nb] = float(intersections.get(tuple(sorted((src_b, nb_b))), 0.0))
        out[src] = refreshed
    return out

# find pair weight from a refilled block, checking both directions
def find_pair_weight(block: dict, a: str, b: str):
    for src, nb in ((a, b), (b, a)):
        d = block.get(src)
        if isinstance(d, dict) and nb in d:
            try:
                return float(d[nb])
            except Exception:
                return None
    return None

def _fmt(x):
    return "NA" if x is None else f"{x:.6g}"

def main():
    # Load fixed sample topology once
    onehop_sample = json.loads(SAMPLE_ONE_HOP_JSON.read_text())
    twohop_sample = json.loads(SAMPLE_TWO_HOP_JSON.read_text())
    all_sample    = json.loads(SAMPLE_ALL_JSON.read_text())

    needed_pairs = collect_needed_pairs_from_samples(
        onehop_sample,
        twohop_sample,
        all_sample
    )

    inter_files = sorted(INTERSECT_DIR.glob("ko_intersection_*.tsv"))
    if not inter_files:
        raise FileNotFoundError(f"No files in {INTERSECT_DIR}")

    for inter_path in inter_files:
        full_base = inter_path.stem
        short_base = full_base.replace("ko_intersection_ko_matrix_sampleids_", "")
        counts_name = full_base.replace("ko_intersection_", "ko_counts_") + ".tsv"
        counts_path = COUNTS_DIR / counts_name

        print(f"\n=== Refill: {short_base} ===")
        if not counts_path.exists():
            print(f"WARNING: counts not found: {counts_path} (treat missing counts as 0)")

        counts = load_counts(counts_path)
        intersections = load_intersections_sparse(inter_path, needed_pairs)

        # Refill all three
        one_refilled = refill_weights(refill_counts(onehop_sample, counts), intersections)
        two_refilled = refill_weights(refill_counts(twohop_sample, counts), intersections)
        all_refilled = refill_weights(refill_counts(all_sample,  counts), intersections)

        # --- DEBUG PRINTS ---
        mv_00150 = get_matrix_value(inter_path, "K00150", "K15634")
        mv_11389 = get_matrix_value(inter_path, "K11389", "K15634")

        oh_00150 = find_pair_weight(one_refilled, "K00150", "K15634")
        oh_11389 = find_pair_weight(one_refilled, "K11389", "K15634")

        th_00150 = find_pair_weight(two_refilled, "K00150", "K15634")
        th_11389 = find_pair_weight(two_refilled, "K11389", "K15634")

        al_00150 = find_pair_weight(all_refilled, "K00150", "K15634")
        al_11389 = find_pair_weight(all_refilled, "K11389", "K15634")

        print(f"[{short_base}] Matrix  row=K00150 col=K15634 -> {_fmt(mv_00150)}")
        print(f"[{short_base}] Matrix  row=K11389 col=K15634 -> {_fmt(mv_11389)}")
        print(f"[{short_base}] OneHop  K00150—K15634 -> {_fmt(oh_00150)} ; K11389—K15634 -> {_fmt(oh_11389)}")
        print(f"[{short_base}] TwoHop  K00150—K15634 -> {_fmt(th_00150)} ; K11389—K15634 -> {_fmt(th_11389)}")
        print(f"[{short_base}] AllHop  K00150—K15634 -> {_fmt(al_00150)} ; K11389—K15634 -> {_fmt(al_11389)}")

        # --- WRITE OUTPUTS ---
        out_one = OUT_ONE_DIR / f"One_Hop_Refilled_{short_base}.json"
        out_two = OUT_TWO_DIR / f"Two_Hop_Refilled_{short_base}.json"
        out_all = OUT_ALL_DIR / f"Module_AllHop_Refilled_{short_base}.json"

        out_one.write_text(json.dumps(one_refilled, indent=2))
        out_two.write_text(json.dumps(two_refilled, indent=2))
        out_all.write_text(json.dumps(all_refilled, indent=2))

        print(f"Wrote:\n  {out_one}\n  {out_two}\n  {out_all}")

    print("\nDone: all lineage-specific outputs generated.")


if __name__ == "__main__":
    main()
