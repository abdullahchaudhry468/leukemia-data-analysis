from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

project_dir = Path(__file__).resolve().parent
data_dir = project_dir / "data"
figures_dir = project_dir / "figures"

figures_dir.mkdir(exist_ok=True)

expression = pd.read_csv(
    data_dir / "expression_normalized.csv",
    index_col=0
)

metadata = pd.read_csv(
    data_dir / "sample_metadata.csv",
    index_col=0
)

counts = metadata["ALL.AML"].value_counts()

plt.figure(figsize=(7, 5))
counts.plot(
    kind="bar",
    color=["steelblue", "indianred"]
)

plt.title("Leukemia Sample Distribution")
plt.xlabel("Leukemia Type")
plt.ylabel("Number of Samples")
plt.xticks(rotation=0)

for i, count in enumerate(counts):
    plt.text(
        i,
        count + 0.5,
        str(count),
        ha="center"
    )

plt.tight_layout()
plt.savefig(
    figures_dir / "sample_distribution.png",
    dpi=300
)
plt.close()

X = expression.T.copy()
X = X.replace([np.inf, -np.inf], np.nan)
X = X.fillna(X.median()).fillna(0)
X = X.loc[:, X.var(axis=0) > 0]

sample_lookup = metadata.copy()
sample_lookup["Samples"] = sample_lookup["Samples"].astype(str)
sample_lookup = sample_lookup.set_index("Samples")

labels = sample_lookup.loc[
    X.index,
    "ALL.AML"
].to_numpy()

X_scaled = StandardScaler().fit_transform(X)

pca = PCA(n_components=2)
components = pca.fit_transform(X_scaled)

plt.figure(figsize=(8, 6))

for leukemia_type, color in [
    ("ALL", "steelblue"),
    ("AML", "indianred")
]:
    mask = labels == leukemia_type

    plt.scatter(
        components[mask, 0],
        components[mask, 1],
        label=leukemia_type,
        color=color,
        alpha=0.8,
        s=65
    )

plt.xlabel(
    f"PC1 ({pca.explained_variance_ratio_[0] * 100:.1f}% variance)"
)

plt.ylabel(
    f"PC2 ({pca.explained_variance_ratio_[1] * 100:.1f}% variance)"
)

plt.title("PCA of Golub Leukemia Gene Expression")
plt.legend(title="Leukemia type")
plt.tight_layout()

plt.savefig(
    figures_dir / "pca_leukemia.png",
    dpi=300
)

plt.close()

print(
    "PCA variance explained:",
    pca.explained_variance_ratio_
)
