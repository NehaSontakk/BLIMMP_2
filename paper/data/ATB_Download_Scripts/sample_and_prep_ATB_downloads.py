import pandas as pd
import os

# High quality sample list from OSF
hq_df = pd.read_csv(
    "hq_seq.samples.txt",
    header=None,
    dtype=str
)
hq_ids = hq_df[0].tolist()

# Full file list with locations
df = pd.read_csv(
    "file_list.all.20240805.tsv",
    sep="\t",
    dtype=str,
    comment="#",
    na_filter=False
)

# Get only high quality samples
sample_col = df.columns[0]
df = df[df[sample_col].isin(hq_ids)]

# 3. Stratified sampling up to 5 per species_sylph
filtered = (
    df
    .groupby("species_sylph",group_keys=False)
    .apply(lambda grp: grp.sample(n=min(len(grp), 5), random_state=42))
    .reset_index(drop=True)
)

filtered.to_csv(
    "filtered_samples_filelist_080725.tsv",
    sep="\t",
    index=False
)

os.makedirs("ATB_Downloads", exist_ok=True)

# wget commands to txt
with open("download_commands.txt", "w") as out_wget:
    for _, row in filtered.iterrows():
        tar = row["tar_xz"]
        url = row["tar_xz_url"]
        out_wget.write("wget -O ATB_Downloads/{0} {1}\n".format(tar, url))

# unzip commands to txt
with open("unzip_commands.txt", "w") as out_unzip:
    for _, row in filtered.iterrows():
        tar = row["tar_xz"]
        out_unzip.write("tar -xJf ATB_Downloads/{0}\n".format(tar))

print("Done")

