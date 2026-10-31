# BEN 2026-11-01
# Change aim: Cluster PBMC single-cell RNA-seq profiles and find marker genes per cluster.
# File: 01_cluster_and_markers.py
"""Step 1: PCA -> neighbour graph -> UMAP -> Leiden clustering -> marker genes.
Data: scanpy.datasets.pbmc68k_reduced() (700 cells x 765 genes, 10x PBMC, bundled with Scanpy)."""
import os
import warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

os.makedirs("results", exist_ok=True)
SEED = 0

# BEN 2026-11-01
# Change aim: Load the bundled PBMC subset with log-normalised expression as the main matrix.
# File: 01_cluster_and_markers.py
# Function: load_pbmc
def load_pbmc():
    a = sc.datasets.pbmc68k_reduced()
    if a.raw is None:
        raise RuntimeError("This Scanpy preset has no .raw; use a default Scanpy preset (1.10 to 1.12).")
    adata = a.raw.to_adata()              # log-normalised expression (not scaled)
    adata.obs = a.obs.copy()
    return adata

adata = load_pbmc()
print(adata)
adata.layers["lognorm"] = adata.X.copy()

# --- quality overview (values are precomputed in the dataset)
fig, ax = plt.subplots(1, 3, figsize=(10, 3))
for a_, k in zip(ax, ["n_genes", "n_counts", "percent_mito"]):
    a_.hist(adata.obs[k], bins=40)
    a_.set_title(k)
plt.tight_layout(); plt.savefig("results/qc_distributions.png", dpi=150); plt.close()

# --- scale, PCA, neighbours, UMAP, Leiden
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, n_comps=30, random_state=SEED)
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=20, random_state=SEED)
sc.tl.umap(adata, random_state=SEED)
sc.tl.leiden(adata, resolution=0.8, random_state=SEED, flavor="igraph", n_iterations=2, directed=False)
print("clusters:", adata.obs["leiden"].nunique())

# --- how well do unsupervised clusters match the reference labels?
ari = adjusted_rand_score(adata.obs["bulk_labels"], adata.obs["leiden"])
nmi = normalized_mutual_info_score(adata.obs["bulk_labels"], adata.obs["leiden"])
print(f"ARI={ari:.3f}  NMI={nmi:.3f}")
with open("results/clustering_summary.txt", "w") as f:
    f.write(f"cells={adata.n_obs} genes={adata.n_vars}\n")
    f.write(f"leiden_clusters={adata.obs['leiden'].nunique()}  resolution=0.8\n")
    f.write(f"ARI_vs_bulk_labels={ari:.3f}\nNMI_vs_bulk_labels={nmi:.3f}\n")

ct = pd.crosstab(adata.obs["leiden"], adata.obs["bulk_labels"])
ct.to_csv("results/cluster_vs_label_table.csv")

sc.pl.umap(adata, color=["leiden", "bulk_labels"], show=False, wspace=0.6)
plt.savefig("results/umap_leiden_and_labels.png", dpi=150, bbox_inches="tight"); plt.close()

# --- marker genes per cluster (Wilcoxon rank-sum on log-normalised expression)
adata.X = adata.layers["lognorm"]
sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon", use_raw=False)
top = sc.get.rank_genes_groups_df(adata, group=None)
top = top.groupby("group", observed=True).head(5)
top.to_csv("results/top5_markers_per_cluster.csv", index=False)

known = ["CD3D", "CD3E", "IL7R", "CD8A", "NKG7", "GNLY", "CD79A", "MS4A1", "CD14", "LYZ",
         "CST3", "FCER1A", "FCGR3A", "S100A8", "S100A9", "PPBP"]
present = [g for g in known if g in adata.var_names]
print("known markers present:", present)
dp = sc.pl.dotplot(adata, present, groupby="leiden", return_fig=True, use_raw=False)
dp.savefig("results/dotplot_known_markers.png", dpi=150)
