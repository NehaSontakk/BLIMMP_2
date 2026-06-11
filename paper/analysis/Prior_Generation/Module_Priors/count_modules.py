import pandas as pd

df = pd.read_csv(
    "metabolic_sample_module_lastcol.tsv",
    sep="\t",
    usecols=["Sample_ID", "module_id", "presence"]
)

df = df[df["module_id"].str.match(r"^M\d{5}$")]
df = df[df["presence"].isin(["Present", "Absent"])]

df = df.drop_duplicates(subset=["Sample_ID", "module_id"])

n_genomes = df["Sample_ID"].nunique()
print(f"n_genomes = {n_genomes}")

summary = (
    df.groupby("module_id")
      .agg(
          present=("presence", lambda x: (x == "Present").sum()),
          total=("presence", "size"),
      )
)
summary["fraction"] = summary["present"] / n_genomes

summary[["present", "fraction"]].rename(
    columns={"fraction": "fraction"}
).to_csv("module_presence_fraction.tsv", sep="\t", index=True)

with open("module_freq.txt", "w") as f:
    for module_id, row in summary.iterrows():
        f.write(f"{module_id}\t{row['fraction']}\n")

print(f"n_modules = {len(summary)}")
print(f"Done. module_freq.txt and module_presence_fraction.tsv written.")
