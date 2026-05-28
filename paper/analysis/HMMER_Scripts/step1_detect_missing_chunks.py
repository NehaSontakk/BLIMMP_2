import os
import glob

def main():
    # Path to the completed HMMER output files
    completed_files = glob.glob('/xdisk/cgoubert/nsontakke/ATB_HMMER_Run2/*.domtblout')
    completed_pairs = set()
    for f in completed_files:
        base = os.path.basename(f)
        name, _ = os.path.splitext(base)
        # Expect format: SAMPLE_ORFs_chunkXX
        sample_part, chunk_part = name.split('_chunk')
        sample = sample_part.replace('_ORFs', '')
        try:
            chunk = int(chunk_part)
        except ValueError:
            continue
        completed_pairs.add((sample, chunk))

    # Get list of all genomes
    genome_files = glob.glob('/xdisk/twheeler/nsontakke/ATB_Analysis_0725/ATB_PRODIGAL/*.fna')
    genomes = [os.path.basename(f).split('_ORFs.fna')[0] for f in genome_files]

    total_genomes = len(genomes)
    completed_genomes = len({sample for sample, _ in completed_pairs})
    completed_chunks = len(completed_pairs)
    expected_chunks = total_genomes * 130

    print(f'Total genomes: {total_genomes}')
    print(f'Genomes with at least one completed chunk: {completed_genomes}')
    print(f'Completed chunks: {completed_chunks} / {expected_chunks}')

    # Generate list of missing chunks
    missing = []
    for sample in genomes:
        for chunk in range(1, 131):
            if (sample, chunk) not in completed_pairs:
                missing.append(f'{sample}_ORFs_chunk{chunk}.domtblout')

    # Write missing list to file
    with open('missing_chunks.txt', 'w') as out:
        for item in missing:
            out.write(f'{item}\n')

    print(f'\nMissing chunks written to missing_chunks.txt ({len(missing)} entries)')

if __name__ == '__main__':
    main()
