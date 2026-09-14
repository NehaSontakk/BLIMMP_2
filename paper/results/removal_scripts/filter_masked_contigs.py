#!/usr/bin/env python3
"""
Drop contigs that are fully or near-fully N-masked from a spliced genome
before handing it to DRAM.

DRAM's rRNA step (barrnap -> nhmmer) fails outright when it hits a contig
with no real nucleotide signal to detect an alphabet from:
    "Invalid alphabet type in target for nhmmer. Expect DNA or RNA."
This is about composition, not length, so DRAM's own --min_contig_size
flag doesn't help; a short-but-real contig is fine, a long-but-all-N
contig is not. This script removes contigs whose N fraction meets or
exceeds a threshold, and writes a report of what got dropped and why.

Usage:
    python filter_masked_contigs.py \
        -i /path/to/GENOME_spliced.fna \
        -o /path/to/GENOME_filtered.fna \
        -r /path/to/GENOME_filtered_report.tsv \
        -t 0.99
"""

import argparse


def read_fasta(filepath):
    """Parse a FASTA file into an ordered list of (header, sequence) tuples."""
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


def filter_contigs(records, threshold):
    """
    Split records into (kept, dropped), where dropped contigs have an N
    fraction >= threshold. Each dropped entry carries the stats needed for
    the report: (header, length, n_count, n_fraction).
    """
    kept = []
    dropped = []
    for header, seq in records:
        length = len(seq)
        if length == 0:
            dropped.append((header, length, 0, 1.0))
            continue
        n_count = seq.upper().count("N")
        n_frac = n_count / length
        if n_frac >= threshold:
            dropped.append((header, length, n_count, n_frac))
        else:
            kept.append((header, seq))
    return kept, dropped


def run_filter(input_path, output_path, report_path, threshold):
    records = read_fasta(input_path)
    if not records:
        raise ValueError(f"No sequences found in {input_path}")

    kept, dropped = filter_contigs(records, threshold)

    if not kept:
        raise ValueError(
            f"All {len(records)} contig(s) in {input_path} met or exceeded the "
            f"N-fraction threshold ({threshold}); nothing left to write. Check "
            f"the input, or lower the threshold if this is unexpected."
        )

    write_fasta(kept, output_path)

    with open(report_path, "w") as f:
        f.write("contig_id\tlength\tn_count\tn_fraction\tstatus\n")
        for header, length, n_count, n_frac in dropped:
            f.write(f"{header}\t{length}\t{n_count}\t{n_frac:.4f}\tdropped\n")

    total = len(records)
    print(
        f"Kept {len(kept)}/{total} contig(s), dropped {len(dropped)} "
        f"(N fraction >= {threshold}) -> {output_path}"
    )


def parse_args():
    parser = argparse.ArgumentParser(
        description="Drop fully or near-fully N-masked contigs from a spliced "
                    "genome so DRAM's rRNA step doesn't choke on them."
    )
    parser.add_argument("-i", "--input", required=True, help="Spliced genome .fna")
    parser.add_argument("-o", "--output", required=True, help="Filtered .fna to write")
    parser.add_argument("-r", "--report", required=True, help="TSV report of dropped contigs")
    parser.add_argument(
        "-t", "--n-fraction-threshold",
        type=float,
        default=0.99,
        help="Contigs with N fraction >= this value are dropped (default 0.99, "
             "i.e. essentially all-N only). Lower this if DRAM still fails on "
             "near-empty contigs after filtering at the default."
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    run_filter(args.input, args.output, args.report, args.n_fraction_threshold)
