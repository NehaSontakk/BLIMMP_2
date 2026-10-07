#!/usr/bin/env python3
"""
Compare old (buggy) vs new (fixed) CoOccurrence_AllVsAll JSON files for a
target set of KOs, extracting each in a single streaming pass per file
(no full-file loading -- both files are multi-GB).

For every changed value, checks whether the difference is an exact
multiple of 256 -- the expected signature of the uint8 overflow bug.
Any changed value that ISN'T a multiple of 256 is flagged as unexpected
and worth investigating separately.

Usage:
    python compare_cooccurrence_old_vs_new.py \
        --old OLD_FILE.json --new NEW_FILE.json \
        --kos K02519 K00001 K00002 ...
"""
import argparse
import sys


def extract_target_blocks(path, target_kos):
    """Single pass through the file, collecting every target KO's block."""
    targets_remaining = set(target_kos)
    results = {}
    current_ko = None
    current_block = None

    with open(path) as f:
        for line in f:
            stripped = line.strip()
            if stripped.endswith(': {'):
                key = stripped.split('"')[1]
                if key in targets_remaining:
                    current_ko = key
                    current_block = {}
                else:
                    current_ko = None
                continue
            if current_ko is None:
                continue
            if stripped in ('},', '}'):
                results[current_ko] = current_block
                targets_remaining.discard(current_ko)
                current_ko = None
                if not targets_remaining:
                    break
                continue
            try:
                key_part, val_part = stripped.split(':', 1)
                other_ko = key_part.strip().strip('"')
                val = float(val_part.strip().rstrip(','))
                current_block[other_ko] = val
            except ValueError:
                continue

    return results


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--old", required=True)
    ap.add_argument("--new", required=True)
    ap.add_argument("--kos", nargs="+", required=True)
    args = ap.parse_args()

    print(f"Extracting {len(args.kos)} KO blocks from OLD file...", file=sys.stderr)
    old_blocks = extract_target_blocks(args.old, args.kos)
    print(f"Extracting {len(args.kos)} KO blocks from NEW file...", file=sys.stderr)
    new_blocks = extract_target_blocks(args.new, args.kos)

    total_compared = 0
    total_changed = 0
    total_unexpected = 0
    examples_shown = 0

    for ko in args.kos:
        old_entry = old_blocks.get(ko)
        new_entry = new_blocks.get(ko)
        if old_entry is None or new_entry is None:
            print(f"{ko}: missing from {'OLD' if old_entry is None else 'NEW'} file, skipping")
            continue

        all_neighbors = set(old_entry) | set(new_entry)
        n_changed_here = 0
        for nb in sorted(all_neighbors):
            old_v = old_entry.get(nb)
            new_v = new_entry.get(nb)
            total_compared += 1
            if old_v != new_v:
                total_changed += 1
                n_changed_here += 1
                diff = (new_v - old_v) if (old_v is not None and new_v is not None) else None
                is_mult_256 = (diff is not None and diff > 0 and diff % 256 == 0)
                if not is_mult_256:
                    total_unexpected += 1
                if examples_shown < 15:
                    tag = "" if is_mult_256 else "  <-- NOT a multiple of 256, unexpected!"
                    print(f"{ko} x {nb}: old={old_v} new={new_v} diff={diff}{tag}")
                    examples_shown += 1

        print(f"{ko}: {n_changed_here}/{len(all_neighbors)} neighbor values changed")

    print(f"\n=== Summary ===")
    print(f"Total value-pairs compared: {total_compared}")
    print(f"Total changed: {total_changed}")
    print(f"Changed values that are NOT a clean multiple of 256: {total_unexpected}")
    if total_unexpected == 0 and total_changed > 0:
        print("\nAll changes are clean multiples of 256 -- fully consistent with the "
              "uint8 overflow bug as the sole explanation for every difference found.")
    elif total_unexpected > 0:
        print("\nSome changes do NOT fit the simple mod-256 pattern -- worth looking "
              "into these specific cases further before assuming the fix is complete.")


if __name__ == "__main__":
    main()
