"""One-time export: adata_qc.h5ad -> Matrix Market + CSVs for the R/Seurat twin.
Run in the endotypes env. Output goes to the data folder (outside git)."""
import os, gzip
import pandas as pd
import scanpy as sc
import scipy.io

src = os.path.expanduser("~/Documents/immune-endotypes/data/smillie_uc/adata_qc.h5ad")
adata = sc.read_h5ad(src)
adata.obs = adata.obs.rename(columns={"Health": "condition", "Location": "site",
                                      "Cluster": "author_label"})
out = os.path.expanduser("~/Documents/immune-endotypes/data/smillie_uc/r_export")
os.makedirs(out, exist_ok=True)
with gzip.open(os.path.join(out, "counts_genes_by_cells.mtx.gz"), "wb") as fh:
    scipy.io.mmwrite(fh, adata.layers["counts"].T.tocsc())
adata.obs.to_csv(os.path.join(out, "obs.csv"))
pd.DataFrame({"gene": adata.var_names}).to_csv(os.path.join(out, "genes.csv"), index=False)
pd.DataFrame({"barcode": adata.obs_names}).to_csv(os.path.join(out, "barcodes.csv"), index=False)
print("exported", adata.shape, "->", out)
