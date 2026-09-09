#!/usr/bin/env python3
"""
Splice out (mask) the genomic regions corresponding to previously-removed
ORFs, so the resulting .fna resembles an incomplete genome recovered from
the environment. Masked with N's.

Usage:
    python splice_orfs_from_genome.py \
        -f /path/to/GENOME.fna \
        -a /path/to/GENOME_ORFs.faa \
        -s /path/to/GENOME_SUBSAMPLES \
        -n GENOME_NAME \
        -o /path/to/OUTPUT_DIR
"""

import os
import glob
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


def parse_orf_coordinates(orfs_faa_path):
    coords = {}
    with open(orfs_faa_path, "r") as f:
        for line in f:
            if not line.startswith(">"):
                continue
            line = line.rstrip()
            fields = line[1:].split(" # ")
            if len(fields) < 3:
                raise ValueError(f"Unexpected Prodigal header format: {line}")
            orf_id = fields[0].strip()
            start = int(fields[1].strip())
            end = int(fields[2].strip())
            contig_id, _, gene_num = orf_id.rpartition("_")
            if not contig_id or not gene_num.isdigit():
                raise ValueError(
                    f"Could not split ORF_ID into contig and gene number: {orf_id}"
                )
            coords[orf_id] = (contig_id, start, end)
    return coords


def find_removed_csvs(subsamples_dir, genome_name):
    """
    Locate every removed-ORF CSV produced by subsample_orfs.py for this
    genome, one per (percent, replicate) combination.
    """
    pattern = os.path.join(
        subsamples_dir, f"{genome_name}_*percentremoved_replicate*.csv"
    )
    return sorted(glob.glob(pattern))


def parse_removed_ids(csv_path):
    with open(csv_path, "r") as f:
        lines = [line.strip() for line in f.readlines()[1:] if line.strip()]
    return lines


def merge_intervals(intervals):
    """Merge overlapping or adjacent (start, end) 1-based inclusive intervals."""
    if not intervals:
        return []
    intervals = sorted(intervals)
    merged = [intervals[0]]
    for start, end in intervals[1:]:
        last_start, last_end = merged[-1]
        if start <= last_end + 1:
            merged[-1] = (last_start, max(last_end, end))
        else:
            merged.append((start, end))
    return merged


def splice_genome(genome_records, contig_intervals, tag):
    """
    Return a new list of (header, sequence) records with the given
    per-contig intervals masked to N. Also returns the count of bases
    masked, for verification against the merged interval lengths.
    """
    spliced = []
    total_masked = 0
    for header, seq in genome_records:
        contig_id = header.split()[0]
        intervals = contig_intervals.get(contig_id, [])
        if not intervals:
            spliced.append((header, seq))
            continue
        seq_chars = list(seq)
        for start, end in intervals:
            if start < 1 or end > len(seq_chars) or start > end:
                raise ValueError(
                    f"[{tag}] Interval ({start},{end}) out of bounds for "
                    f"contig {contig_id} (length {len(seq_chars)})"
                )
            span = end - start + 1
            seq_chars[start - 1:end] = ["N"] * span
            total_masked += span
        spliced.append((header, "".join(seq_chars)))
    return spliced, total_masked


def verify_masked_count(total_masked, contig_intervals, expected_orf_bases, tag):
    """
    Sanity check: the number of bases masked should equal the sum of the
    merged interval lengths, not the raw sum of removed ORF lengths, since a
    small number of Prodigal ORFs can overlap at the edges and merging
    collapses that overlap. Report both so a large mismatch is easy to spot.
    """
    merged_bases = sum(
        end - start + 1 for ivs in contig_intervals.values() for start, end in ivs
    )
    if total_masked != merged_bases:
        raise AssertionError(
            f"[{tag}] Masked base count ({total_masked}) does not match "
            f"merged interval bases ({merged_bases})"
        )
    if expected_orf_bases and total_masked > expected_orf_bases:
        raise AssertionError(
            f"[{tag}] Masked more bases ({total_masked}) than the raw sum of "
            f"removed ORF lengths ({expected_orf_bases}). Merging intervals "
            f"should only ever reduce this number, never increase it."
        )


def splice_all(fna_path, orfs_faa_path, subsamples_dir, genome_name, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    genome_records = read_fasta(fna_path)
    orf_coords = parse_orf_coordinates(orfs_faa_path)

    csv_paths = find_removed_csvs(subsamples_dir, genome_name)
    if not csv_paths:
        raise FileNotFoundError(
            f"No removed-ORF CSVs found for {genome_name} in {subsamples_dir}. "
            f"Expected pattern {genome_name}_<pct>percentremoved_replicate<r>.csv"
        )

    print(f"[{genome_name}] Loaded genome with {len(genome_records)} contig(s) from {fna_path}")
    print(f"[{genome_name}] Loaded {len(orf_coords)} ORF coordinate(s) from {orfs_faa_path}")
    print(f"[{genome_name}] Found {len(csv_paths)} removed-ORF CSV(s) to splice")

    for csv_path in csv_paths:
        tag = os.path.splitext(os.path.basename(csv_path))[0]
        out_path = os.path.join(output_dir, f"{tag}_spliced.fna")

        if os.path.exists(out_path):
            print(f"  [{genome_name}] SKIP {tag} (already spliced)")
            continue

        removed_ids = parse_removed_ids(csv_path)

        contig_raw_intervals = {}
        expected_orf_bases = 0
        missing_ids = []
        for orf_id in removed_ids:
            if orf_id not in orf_coords:
                missing_ids.append(orf_id)
                continue
            contig_id, start, end = orf_coords[orf_id]
            contig_raw_intervals.setdefault(contig_id, []).append((start, end))
            expected_orf_bases += end - start + 1

        if missing_ids:
            raise ValueError(
                f"[{tag}] {len(missing_ids)} removed ORF_ID(s) from the CSV were not "
                f"found in {orfs_faa_path}, e.g. {missing_ids[:5]}. The .faa used for "
                f"coordinates must be the exact same one subsample_orfs.py drew from."
            )

        contig_intervals = {
            contig_id: merge_intervals(ivs) for contig_id, ivs in contig_raw_intervals.items()
        }

        spliced_records, total_masked = splice_genome(genome_records, contig_intervals, tag)
        verify_masked_count(total_masked, contig_intervals, expected_orf_bases, tag)

        write_fasta(spliced_records, out_path)

        genome_len = sum(len(seq) for _, seq in genome_records)
        pct_masked = 100.0 * total_masked / genome_len if genome_len else 0.0
        print(
            f"  [{genome_name}] {tag}: masked {total_masked} bases across "
            f"{len(contig_intervals)} contig(s) ({pct_masked:.2f}% of genome) "
            f"-> {os.path.basename(out_path)} [verified OK]"
        )


def parse_args():
    parser = argparse.ArgumentParser(
        description="Splice (N-mask) the genomic regions of previously-removed "
                    "ORFs to simulate an incomplete recovered genome."
    )
    parser.add_argument("-f", "--fna", required=True, help="Original genome .fna")
    parser.add_argument("-a", "--orfs-faa", required=True, help="Prodigal .faa used to select ORFs for removal")
    parser.add_argument("-s", "--subsamples-dir", required=True, help="Dir containing removed-ORF CSVs from subsample_orfs.py")
    parser.add_argument("-n", "--name", required=True, help="Genome name, must match the CSV filename prefix")
    parser.add_argument("-o", "--output-dir", required=True, help="Directory to write spliced .fna files")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    splice_all(args.fna, args.orfs_faa, args.subsamples_dir, args.name, args.output_dir)
