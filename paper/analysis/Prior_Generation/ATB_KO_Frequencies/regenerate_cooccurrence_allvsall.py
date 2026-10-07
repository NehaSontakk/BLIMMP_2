#!/usr/bin/env python3
"""
Regenerate CoOccurrence_AllVsAll_domain_level_priors.json from the current
ko_intersection_ko_matrix_sampleids_domain_level_priors.tsv, which is now
internally consistent with the rebuilt ko_matrix.tsv/counts (post K99001,
post the changed-calls adjudication fixes).

The OLD CoOccurrence_AllVsAll file predates the matrix rebuild entirely, and
is now inconsistent with current counts for KOs whose presence calls changed
during adjudication -- not just K99001. That inconsistency causes a hard
crash in BLIMMP's SubstitutionPipeline (Nij > Nj check). This fully replaces
it with data drawn from the same source as everything else currently in
production.

Streams the input TSV line-by-line and writes the output JSON line-by-line
(hand-rolled, not json.dump on the full structure) to stay memory-safe
regardless of the matrix's size -- same philosophy as
insert_mhpp001_into_module_skeleton.py's streaming submatrix loader, just
applied to the full matrix instead of a 9x9 slice.

Only non-zero entries are written per KO (matching the sparse nature of the
source data and keeping the output file size reasonable), plus a "_count"
field using diagonal self-intersection (a KO's intersection with itself
equals its own prevalence count).

Usage:
    python regenerate_cooccurrence_allvsall.py \
        --intersection-tsv .../ko_intersection_ko_matrix_sampleids_domain_level_priors.tsv \
        --out .../CoOccurrence_AllVsAll_domain_level_priors.json
"""
import argparse
import json
import sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--intersection-tsv", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--progress-every", type=int, default=2000)
    args = ap.parse_args()

    n_kos_written = 0
    n_edges_written = 0

    with open(args.intersection_tsv) as fin, open(args.out, "w") as fout:
        header = fin.readline().rstrip("\n").split("\t")
        col_names = header[1:]  # KO column labels, in order
        n_cols = len(col_names)
        col_name_to_idx = {name: i for i, name in enumerate(col_names)}
        print(f"KO universe size: {n_cols}", file=sys.stderr)

        fout.write("{\n")
        first_ko = True

        for lineno, line in enumerate(fin, 1):
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            row_ko = parts[0]
            values = parts[1:]

            if len(values) != n_cols:
                print(f"[WARN] line {lineno} ({row_ko}): {len(values)} values, "
                      f"expected {n_cols} -- skipping", file=sys.stderr)
                continue

            # self-intersection = this KO's own prevalence count
            self_idx = col_name_to_idx.get(row_ko)
            own_count = float(values[self_idx]) if self_idx is not None else None

            entry_parts = []
            for col_name, v in zip(col_names, values):
                if col_name == row_ko:
                    continue
                if v and v != "0" and v != "0.0":
                    fv = float(v)
                    if fv != 0.0:
                        entry_parts.append(f'    "{col_name}": {fv}')
                        n_edges_written += 1

            if not first_ko:
                fout.write(",\n")
            first_ko = False

            fout.write(f'  "{row_ko}": {{\n')
            if own_count is not None:
                fout.write(f'    "_count": {own_count}')
                if entry_parts:
                    fout.write(",\n")
            if entry_parts:
                fout.write(",\n".join(entry_parts))
            fout.write("\n  }")

            n_kos_written += 1
            if n_kos_written % args.progress_every == 0:
                print(f"  ...{n_kos_written} KOs written", file=sys.stderr)

        fout.write("\n}\n")

    print(f"\nDone. Wrote {n_kos_written} KOs, {n_edges_written} non-zero edges "
          f"-> {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
