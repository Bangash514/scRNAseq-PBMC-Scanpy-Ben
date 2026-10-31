# BEN 2026-11-01
# Change aim: Train and cross-validate a cell-type classifier on PBMC gene expression.
# File: 02_celltype_classifier.py
"""Step 2: supervised cell-type classification (logistic regression vs random forest)
with stratified 5-fold cross-validation, confusion matrix and top genes per class."""
import os
import warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, f1_score, confusion_matrix, ConfusionMatrixDisplay
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

os.makedirs("results", exist_ok=True)
SEED = 0

a = sc.datasets.pbmc68k_reduced()
ad = a.raw.to_adata()
X = np.asarray(ad.X.todense()) if hasattr(ad.X, "todense") else np.asarray(ad.X)
y = a.obs["bulk_labels"].astype(str).values
genes = np.array(ad.var_names)
classes = sorted(set(y))
print("cells, genes:", X.shape, " classes:", len(classes))

cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=SEED)
models = {
    "Logistic regression": make_pipeline(StandardScaler(), LogisticRegression(
        C=0.01, max_iter=3000, class_weight="balanced")),
    "Random forest": RandomForestClassifier(n_estimators=300, class_weight="balanced",
                                            n_jobs=-1, random_state=SEED),
}
rows = []
preds = {}
major = pd.Series(y).value_counts(normalize=True).iloc[0]
rows.append(("Majority-class baseline", major, np.nan))
for name, m in models.items():
    p = cross_val_predict(m, X, y, cv=cv)
    preds[name] = p
    rows.append((name, accuracy_score(y, p), f1_score(y, p, average="macro")))
res = pd.DataFrame(rows, columns=["model", "cv_accuracy", "cv_macro_f1"])
print(res)
res.to_csv("results/classifier_cv_metrics.csv", index=False, float_format="%.3f")

best = res.iloc[1:].sort_values("cv_macro_f1").iloc[-1]["model"]
cm = confusion_matrix(y, preds[best], labels=classes, normalize="true")
fig, ax = plt.subplots(figsize=(10, 8.5))
ConfusionMatrixDisplay(cm, display_labels=classes).plot(ax=ax, xticks_rotation=45,
                                                       colorbar=False, values_format=".2f")
plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
ax.set_title(f"{best}: 5-fold CV, row-normalised")
plt.tight_layout(); plt.savefig("results/confusion_cv.png", dpi=150); plt.close()

# top genes per class from a logistic-regression model fitted on all cells
lr = models["Logistic regression"].fit(X, y)
coef = lr[-1].coef_
out = []
for i, c in enumerate(lr[-1].classes_):
    for g in genes[np.argsort(coef[i])[::-1][:5]]:
        out.append((c, g))
pd.DataFrame(out, columns=["cell_type", "top_gene"]).to_csv("results/top_genes_per_class.csv", index=False)
print("done")
