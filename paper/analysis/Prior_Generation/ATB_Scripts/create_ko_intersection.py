import pandas as pd
import numpy as np
from scipy.sparse import csr_matrix

df = pd.read_csv("/xdisk/cgoubert/nsontakke/ATB_KO_Frequencies/ko_matrix.tsv", sep="\t", index_col=0)
X = csr_matrix(df.to_numpy(dtype=np.uint8, copy=False))
M = (X @ X.T).astype(np.int32)

co = pd.DataFrame.sparse.from_spmatrix(M, index=df.index, columns=df.index)

co.to_csv("/xdisk/cgoubert/nsontakke/ATB_KO_Frequencies/ko_intersection.tsv", sep="\t")
