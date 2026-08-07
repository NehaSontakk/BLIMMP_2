#!/usr/bin/env python3
"""
Build an anvi'o external-gene-calls TSV *and* a matching placeholder
nucleotide FASTA directly from a protein FASTA (.faa), so the protein data
can be loaded into anvi-gen-contigs-database without re-calling genes.

IMPORTANT: anvi-gen-contigs-database's --contigs-fasta is always validated
and stored as NUCLEOTIDE sequence, even when an external-gene-calls file
supplies aa_sequence directly. Passing the raw protein FASTA as
--contigs-fasta fails immediately with:
    "Config Error: ... contains characters that are not any of A, C, T, G, N..."
because the protein alphabet (e.g. M, K, V) isn't valid nucleotide sequence.

The correct pattern (confirmed against anvi'o's own docs) is:
  - --contigs-fasta: a placeholder nucleotide sequence for each gene, of
    length 3 * len(protein), covering exactly the gene's start/stop window.
    Its actual bases don't matter, because...
  - ...the external-gene-calls TSV's aa_sequence column tells anvi'o to use
    the given amino acid sequence directly instead of translating the
    nucleotide sequence, per anvi'o's external-gene-calls docs.

Each protein is still treated as its own single-gene "contig" (start=0,
stop=3*aa_len), matching anvi'o's documented convention for externally
provided, already-called genes.

Usage:
    python make_external_gene_calls.py -i sample.faa \
        -o sample_gene_calls.tsv \
        -f sample_placeholder_contigs.fna

Outputs:
  - TSV columns: gene_callers_id, contig, start, stop, direction, partial,
    call_type, source, version, aa_sequence
  - FASTA: one placeholder nucleotide "contig" per protein, headers matching
    the TSV's contig column exactly, length 3*len(protein), filled with 'N'.
"""

import argparse
import sys


def read_fasta(path):
    recs, h, buf = [], None, []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if h is not None:
                    recs.append((h, "".join(buf)))
                h = line[1:].strip().split()[0]
                buf = []
            else:
                s = line.strip()
                if s:
                    buf.append(s)
        if h is not None:
            recs.append((h, "".join(buf)))
    return recs


def write_placeholder_fasta(records, out_path, line_width=60):
    """Write one placeholder nucleotide 'contig' per protein: length
    3*len(protein), filled with N (a valid nucleotide character anvi'o
    accepts; content is irrelevant since aa_sequence overrides translation)."""
    with open(out_path, "w") as f:
        for header, seq in records:
            nt_len = len(seq) * 3
            placeholder = "N" * nt_len
            f.write(f">{header}\n")
            for i in range(0, len(placeholder), line_width):
                f.write(placeholder[i:i + line_width] + "\n")


def main():
    ap = argparse.ArgumentParser(description="Build anvi'o external-gene-calls TSV + placeholder nucleotide FASTA from a protein FASTA.")
    ap.add_argument("-i", "--faa", required=True, help="Input protein FASTA")
    ap.add_argument("-o", "--output", required=True, help="Output gene-calls TSV path")
    ap.add_argument("-f", "--fasta-output", required=True, help="Output placeholder nucleotide FASTA path")
    args = ap.parse_args()

    records = read_fasta(args.faa)
    if not records:
        print(f"ERROR: no sequences found in {args.faa}", file=sys.stderr)
        sys.exit(1)

    with open(args.output, "w") as out:
        out.write(
            "gene_callers_id\tcontig\tstart\tstop\tdirection\tpartial\t"
            "call_type\tsource\tversion\taa_sequence\n"
        )
        for idx, (header, seq) in enumerate(records):
            aa_len = len(seq)
            # anvi'o expects the "contig" to be exactly the length of the gene
            # when treating one protein as its own single-gene contig:
            # start=0, stop=3*len(aa_sequence) (nucleotide-equivalent length).
            # direction must be 'f' or 'r' (forward/reverse) per anvi'o's
            # external-gene-calls spec; we have no strand info from a protein
            # FASTA alone, so 'f' is used uniformly (this is arbitrary but
            # harmless for gene-content/KO annotation purposes since we are
            # supplying the aa_sequence directly and anvi'o will not re-translate).
            stop = aa_len * 3
            out.write(
                f"{idx}\t{header}\t0\t{stop}\tf\t0\t1\tprodigal\tv2.6.3\t{seq}\n"
            )

    write_placeholder_fasta(records, args.fasta_output)

    print(f"Wrote {len(records)} gene calls to {args.output}")
    print(f"Wrote {len(records)} placeholder nucleotide contigs to {args.fasta_output}")


if __name__ == "__main__":
    main()
