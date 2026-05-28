import os
import re
from collections import defaultdict

def generate_command_files(missing_file='missing_chunks.txt', 
                           ids_dir='/xdisk/twheeler/nsontakke/ATB_Analysis_0725/IDS',
                           genome_dir='/xdisk/twheeler/nsontakke/ATB_Analysis_0725/ATB_PRODIGAL',
                           db_dir='/xdisk/twheeler/nsontakke/ATB_Analysis_0725/chunked_hmmdb',
                           output_dir='/xdisk/cgoubert/nsontakke/ATB_HMMER_Run2'):
    with open(missing_file) as mf:
        lines = [line.strip() for line in mf if line.strip()]
    
    # Group commands by chunk ID
    cmds_by_chunk = defaultdict(list)
    pattern = re.compile(r'(.+?)_ORFs_chunk(\d+)\.domtblout')
    
    for entry in lines:
        m = pattern.match(entry)
        if not m:
            continue
        sample, chunk = m.groups()
        chunk = int(chunk)
        tbl = f"{output_dir}/{sample}_chunk{chunk}.tblout"
        dom = f"{output_dir}/{sample}_chunk{chunk}.domtblout"
        out = f"{output_dir}/{sample}_chunk{chunk}.hmmout"
        cmd = (
            f'hmmsearch --cpu 0 '
            f'--domtblout "{dom}" '
            f'-o "{out}" '
            f'"{db_dir}/chunk_{chunk}.hmm" '
            f'"{genome_dir}/{sample}_ORFs.faa"'
        )
        cmds_by_chunk[chunk].append(cmd)
    
    # Ensure IDS directory exists
    
    output_commands = "/xdisk/twheeler/nsontakke/ATB_Analysis_0725/RUN_HMMSCAN_MISSING_FILES/COMMANDS/"
    os.makedirs(output_commands, exist_ok=True)
    # Write commands to files
    for chunk, cmd_list in cmds_by_chunk.items():
        path = os.path.join(output_commands, f"commands_{chunk}.txt")
        with open(path, 'w') as out:
            out.write('\n'.join(cmd_list))
    print(f"Generated {len(cmds_by_chunk)} command files in '{output_commands}'.")

if __name__ == '__main__':
    generate_command_files()
