#!/usr/bin/env python3
import argparse
import os
import glob

def combine_group(files, out_path, header_lines, footer_lines):
    # Read header/footer from first file
    with open(files[0]) as f:
        lines = f.readlines()
    if len(lines) < header_lines + footer_lines:
        raise SystemExit(f"File {files[0]!r} too short for header/footer extraction")
    header = lines[:header_lines]
    footer = lines[-footer_lines:]

    with open(out_path, 'w') as out:
        # Write header
        out.writelines(header)
        # Write bodies of all files
        for fn in files:
            with open(fn) as f:
                all_lines = f.readlines()
            body = all_lines[header_lines:len(all_lines)-footer_lines]
            out.writelines(body)
        # Write footer
        out.writelines(footer)

def main(indir, header_lines, footer_lines):
    pattern = os.path.join(indir, '*_chunk*.domtblout')
    all_files = sorted(glob.glob(pattern))
    if not all_files:
        raise SystemExit(f"No .domtblout files found in {indir!r}")

    # Group files by sample prefix (before _chunk)
    groups = {}
    for filepath in all_files:
        fname = os.path.basename(filepath)
        prefix = fname.split('_chunk')[0]
        groups.setdefault(prefix, []).append(filepath)

    # Combine each group
    for prefix, files in groups.items():
        files_sorted = sorted(files)
        out_file = os.path.join(indir, f"{prefix}.combined.domtblout")
        print(f"Combining {len(files_sorted)} files for '{prefix}' into '{out_file}'")
        combine_group(files_sorted, out_file, header_lines, footer_lines)
    print("All groups combined.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Combine chunked .domtblout files by sample")
    parser.add_argument("-i", "--indir", required=True,
                        help="Directory containing *_chunk*.domtblout files")
    parser.add_argument("-H", "--header", type=int, default=3,
                        help="Number of header lines from first chunk (default: 3)")
    parser.add_argument("-F", "--footer", type=int, default=10,
                        help="Number of footer lines from first chunk (default: 10)")
    args = parser.parse_args()
    main(args.indir, args.header, args.footer)

# Usage example:
# python combine_domtblout_groups.py -i Endosymbiont-genomes-HMMER -H 3 -F 10
