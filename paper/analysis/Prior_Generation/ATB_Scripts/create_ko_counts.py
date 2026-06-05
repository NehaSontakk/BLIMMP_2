import pandas as pd

df = pd.read_csv("/xdisk/cgoubert/nsontakke/ATB_KO_Frequencies/ko_matrix.tsv", sep="\t", index_col=0)
df = (df > 0).astype(int)
ko_counts = df.sum(axis=1).astype(int)
#prevalence = ko_counts / df.shape[1]

out = pd.DataFrame({"KO_count": ko_counts})
out.index.name = "KOs"
out.sort_values("KO_count", ascending=False).to_csv("/xdisk/cgoubert/nsontakke/ATB_KO_Frequencies/ko_counts.tsv", sep="\t")
print(f"Wrote {len(out)} KO counts to ko_counts.tsv")
