# scRNAseq-PBMC-Scanpy-Ben

**Single-cell RNA-seq analysis of human blood immune cells (PBMC) with Scanpy: clustering, marker genes and cell-type classification**

An intermediate bioinformatics project on real single-cell RNA-seq data. It covers the standard transcriptomics workflow (dimensionality reduction, graph-based clustering, marker-gene discovery) and adds a supervised machine-learning step with honest cross-validation.

> Learning and portfolio project. It uses a small, already-preprocessed teaching subset and is not a research result.

## Data
`scanpy.datasets.pbmc68k_reduced()`: 700 cells x 765 genes, a subset of the 10x Genomics PBMC 68k data (Zheng et al., 2017), bundled with Scanpy (no download needed). It ships with 10 cell-type labels (`bulk_labels`), which I use as an approximate reference. Because the data are already filtered to 765 genes and normalised, this project starts at PCA and does not redo raw read-level quality control.

## Pipeline
| Step | Script | What it does |
|---|---|---|
| 1 | `src/01_cluster_and_markers.py` | Scale, PCA (30 PCs), kNN graph, UMAP, Leiden clustering (resolution 0.8), comparison with reference labels (ARI/NMI), Wilcoxon marker genes per cluster, dot plot of known immune markers |
| 2 | `src/02_celltype_classifier.py` | Predict cell type from gene expression: logistic regression vs random forest, stratified 5-fold cross-validation, majority-class baseline, confusion matrix, top genes per class |

## Quick start
```bash
git clone https://github.com/Bangash514/scRNAseq-PBMC-Scanpy-Ben.git
cd scRNAseq-PBMC-Scanpy-Ben
pip install -r requirements.txt
python src/01_cluster_and_markers.py
python src/02_celltype_classifier.py
```
Run from the project root (Python 3.10 or newer). Outputs go to `results/`. Tested with Python 3.13, Scanpy 1.12.5 and scikit-learn 1.9.1; exact numbers can vary slightly with library versions.

## Results (from my run)
**Unsupervised clustering** (`results/clustering_summary.txt`)
- 9 Leiden clusters. Agreement with the reference labels: ARI = 0.524, NMI = 0.657.
- Cluster markers match known immune biology (`results/dotplot_known_markers.png`, `results/top5_markers_per_cluster.csv`):
  - B cells: CD79A, CD79B, MS4A1
  - NK / cytotoxic cells: NKG7, GNLY, GZMB
  - T cells: CD3D, LCK, IL32
  - Dendritic / monocyte-like cells: FCER1A, LYZ, CST3, HLA-DR genes
- One cluster (markers MZB1, IGJ, TNFRSF17) is consistent with plasma-cell-like B lineage cells. This is an interpretation from marker genes only.
- The two T-cell clusters and the CD8 groups overlap in UMAP space and are not cleanly separated by clustering.

**Supervised classification** (`results/classifier_cv_metrics.csv`, 5-fold stratified CV, 10 classes)
| Model | CV accuracy | CV macro-F1 |
|---|---|---|
| Majority-class baseline | 0.343 | n/a |
| Logistic regression | 0.834 | 0.718 |
| Random forest | 0.807 | 0.653 |

- Easy classes (CD34+, B cells, CD14+ monocytes, dendritic cells) are recalled at 0.83 to 1.00 (`results/confusion_cv.png`).
- T-cell subsets are hard: for example CD4+ memory (19 cells) is recalled at 0.37 and CD4+ naive (8 cells) at 0.50. These are very small classes, so these numbers are unstable.
- Top genes learned per class (`results/top_genes_per_class.csv`) include CD79A/MS4A1 (B), GNLY/GZMB (NK), CD8A/CD8B (CD8 naive) and CST3/LYZ (dendritic), which agree with the marker analysis.

Figures: `umap_leiden_and_labels.png`, `dotplot_known_markers.png`, `confusion_cv.png`, `qc_distributions.png`.

## Limitations
- Small dataset (700 cells). Several classes have fewer than 20 cells.
- The reference labels ship with the dataset and their derivation is not documented in the Scanpy docs, so they are an approximate reference, not ground truth. T-cell subtypes overlap strongly in expression.
- Genes were pre-selected on all cells in the shipped dataset, so cross-validation is slightly optimistic. Hyperparameters were not tuned.
- Clustering depends on resolution and random seed (fixed to 0 here).

## Possible extensions
Run the full pipeline on a complete 10x dataset (QC, doublet removal, batch correction), add differential-expression statistics and pathway enrichment, and try an automated annotation tool.

## References
- Wolf F.A. et al., *SCANPY: large-scale single-cell gene expression data analysis*, Genome Biology, 2018.
- Zheng G.X.Y. et al., *Massively parallel digital transcriptional profiling of single cells*, Nature Communications, 2017.
- Traag V.A. et al., *From Louvain to Leiden: guaranteeing well-connected communities*, Scientific Reports, 2019.

## Author
BEN, Software Engineer. GitHub: [Bangash514](https://github.com/Bangash514)
