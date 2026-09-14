#!/usr/bin/env python3
"""
Randomly remove a percentage of ORFs from a Prodigal .faa file, across
multiple percentages and replicates, and write:
  - a CSV listing which ORF_IDs were removed
  - a FASTA of the remaining (kept) ORFs

Usage:
    python subsample_orfs.py -i /path/to/MED4_ORFs.faa -n MED4 -o /path/to/OUTPUT_DIR
"""

import os
import random
import argparse
import sys


def read_fasta(filepath):
    """Parse a FASTA file into a list of (header, sequence) tuples."""
    header = None
    records = []
    seq_lines = []
    with open(filepath, "r") as f:
        for line in f:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_lines)))
                header = line[1:].split()[0]
                seq_lines = []
            else:
                seq_lines.append(line)
        if header is not None:
            records.append((header, "".join(seq_lines)))
    return records


def write_fasta(records, out_path, line_width=60):
    with open(out_path, "w") as f:
        for header, seq in records:
            f.write(f">{header}\n")
            for i in range(0, len(seq), line_width):
                f.write(seq[i:i + line_width] + "\n")


def verify_consistency(all_ids, kept_records, removed_ids, fasta_out, removed_csv, tag):

    #Cross-check that the written FASTA and CSV agree with each other and with the original record set:

    written_kept = read_fasta(fasta_out)
    written_kept_ids = [h for h, _ in written_kept]

    with open(removed_csv, "r") as f:
        lines = [line.strip() for line in f.readlines()[1:] if line.strip()]
    written_removed_ids = lines

    kept_set = set(written_kept_ids)
    removed_set = set(written_removed_ids)

    overlap = kept_set & removed_set
    assert not overlap, (
        f"[{tag}] {len(overlap)} ID(s) appear in BOTH the kept FASTA and the "
        f"removed CSV, e.g. {sorted(overlap)[:5]}"
    )

    # 3: union covers the original ID set exactly
    union = kept_set | removed_set
    missing = all_ids - union
    extra = union - all_ids
    assert not missing, (
        f"[{tag}] {len(missing)} original ID(s) are missing from both outputs, "
        f"e.g. {sorted(missing)[:5]}"
    )
    assert not extra, (
        f"[{tag}] {len(extra)} ID(s) in the outputs were not in the original file, "
        f"e.g. {sorted(extra)[:5]}"
    )

    assert len(written_kept_ids) == len(kept_set), (
        f"[{tag}] Duplicate ID(s) found within the kept FASTA"
    )
    assert len(written_removed_ids) == len(removed_set), (
        f"[{tag}] Duplicate ID(s) found within the removed CSV"
    )


    assert removed_set == set(removed_ids), (
        f"[{tag}] Removed IDs written to CSV don't match the IDs computed in memory"
    )


def subsample_orfs(filepath, genome_name, output_dir, removal_percents, replicates=3):
    records = read_fasta(filepath)
    total_orfs = len(records)
    if total_orfs == 0:
        raise ValueError(f"No sequences found in {filepath}")

    os.makedirs(output_dir, exist_ok=True)

    print(f"[{genome_name}] Loaded {total_orfs} ORFs from {filepath}")

    all_ids = {header for header, _ in records}
    if len(all_ids) != total_orfs:
        raise ValueError(
            f"[{genome_name}] Input FASTA has duplicate headers "
            f"({total_orfs} records but only {len(all_ids)} unique IDs) — "
            f"cannot reliably verify subsampling. Fix duplicate IDs first."
        )

    for percent in removal_percents:
        keep_frac = 1 - (percent / 100.0)
        keep_count = int(round(total_orfs * keep_frac))

        for rep in range(1, replicates + 1):
            # Distinct, reproducible seed per (genome, percent, replicate)
            seed = hash((genome_name, percent, rep)) & 0xFFFFFFFF
            random.seed(seed)

            if keep_count > 0:
                keep_records = random.sample(records, keep_count)
            else:
                keep_records = []
            keep_ids = {header for header, _ in keep_records}
            removed_records = [(h, s) for h, s in records if h not in keep_ids]
            removed_ids = [h for h, _ in removed_records]

            tag = f"{genome_name}_{percent}percentremoved_replicate{rep}"
            fasta_out = os.path.join(output_dir, f"{tag}.faa")
            removed_csv = os.path.join(output_dir, f"{tag}.csv")

            write_fasta(keep_records, fasta_out)
            with open(removed_csv, "w") as f:
                f.write("ORF_ID\n")
                for rid in removed_ids:
                    f.write(f"{rid}\n")

            tag = f"{genome_name}_{percent}pct_rep{rep}"
            verify_consistency(all_ids, keep_records, removed_ids, fasta_out, removed_csv, tag)

            print(
                f"  [{genome_name}] {percent}% removed, rep {rep}: "
                f"kept {len(keep_records)}, removed {len(removed_ids)} "
                f"-> {os.path.basename(fasta_out)} [verified OK]"
            )


def parse_args():
    parser = argparse.ArgumentParser(
        description="Randomly subsample ORFs from a Prodigal .faa file across "
                    "percentages and replicates."
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Path to the input .faa (Prodigal-predicted protein FASTA)"
    )
    parser.add_argument(
        "-n", "--name",
        required=True,
        help="Genome/sample name used as the output file prefix (e.g. MED4)"
    )
    parser.add_argument(
        "-o", "--output-dir",
        required=True,
        help="Directory in which to write subsamples (flat, one folder per genome)"
    )
    parser.add_argument(
        "-p", "--percents",
        nargs="+",
        type=int,
        default=[0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100],
        help="List of percent values to remove (default: 0,10,...,100)"
    )
    parser.add_argument(
        "-r", "--replicates",
        type=int,
        default=3,
        help="Number of replicate subsamples per percent (default: 3)"
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    subsample_orfs(
        args.input,
        args.name,
        args.output_dir,
        args.percents,
        replicates=args.replicates,
    )
