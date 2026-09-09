from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Find the data folder inside this project
project_dir = Path(__file__).resolve().parent
data_dir = project_dir / "data"

# Load the CSV files exported from R
expression = pd.read_csv(data_dir / "expression_normalized.csv", index_col=0)
metadata = pd.read_csv(data_dir / "sample_metadata.csv", index_col=0)

# Inspect the data
print("Expression matrix:", expression.shape)
print("Sample metadata:", metadata.shape)

print("\nFirst five expression rows:")
print(expression.iloc[:5, :5])

print("\nMetadata columns:")
print(metadata.columns.tolist())

print("\nLeukemia types:")
print(metadata["ALL.AML"].value_counts())

# Create a folder for figures
figures_dir = project_dir / "figures"
figures_dir.mkdir(exist_ok=True)

# Count the leukemia types
counts = metadata["ALL.AML"].value_counts()

# Create a bar chart
plt.figure(figsize=(7, 5))
counts.plot(kind="bar", color=["steelblue", "indianred"])

plt.title("Leukemia Sample Distribution")
plt.xlabel("Leukemia Type")
plt.ylabel("Number of Samples")
plt.xticks(rotation=0)

# Add the sample count above each bar
for i, count in enumerate(counts):
    plt.text(i, count + 0.5, str(count), ha="center")

plt.tight_layout()
plt.savefig(figures_dir / "sample_distribution.png", dpi=300)
plt.show()

print("Figure saved to:", figures_dir / "sample_distribution.png")

# ==========================================
# PCA: Explore leukemia expression patterns
# ==========================================

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

# Transpose so rows are samples and columns are genes
X = expression.T.copy()

# Replace infinite values and handle any missing measurements
X = X.replace([np.inf, -np.inf], np.nan)
X = X.fillna(X.median()).fillna(0)

# Remove genes with no variation across samples
X = X.loc[:, X.var(axis=0) > 0]

# Match sample labels to the expression matrix
sample_lookup = metadata.copy()
sample_lookup["Samples"] = sample_lookup["Samples"].astype(str)
sample_lookup = sample_lookup.set_index("Samples")
labels = sample_lookup.loc[X.index, "ALL.AML"].to_numpy()

# Standardize gene measurements
X_scaled = StandardScaler().fit_transform(X)

# Reduce thousands of features to two principal components
pca = PCA(n_components=2)
components = pca.fit_transform(X_scaled)

# Create the PCA figure
plt.figure(figsize=(8, 6))

for leukemia_type, color in [("ALL", "steelblue"), ("AML", "indianred")]:
    mask = labels == leukemia_type
    plt.scatter(
        components[mask, 0],
        components[mask, 1],
        label=leukemia_type,
        color=color,
        alpha=0.8,
        s=65
    )

plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0] * 100:.1f}% variance)")
plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1] * 100:.1f}% variance)")
plt.title("PCA of Golub Leukemia Gene Expression")
plt.legend(title="Leukemia type")
plt.tight_layout()

# Save the figure
plt.savefig(figures_dir / "pca_leukemia.png", dpi=300)
plt.show()

print("PCA figure saved successfully.")
print("Variance explained:", pca.explained_variance_ratio_)

# Data validation before differential expression
print("\n--- DATA VALIDATION ---")
print("Expression shape:", expression.shape)
print("Metadata shape:", metadata.shape)
print("Missing expression values:", expression.isna().sum().sum())
print("Expression minimum:", expression.min().min())
print("Expression maximum:", expression.max().max())
print("Expression median:", expression.stack().median())
print("Sample IDs match:", set(expression.columns.astype(str)) == set(metadata["Samples"].astype(str)))
print("Leukemia types:")
print(metadata["ALL.AML"].value_counts())

# Inspect expression-value distribution
print("\n--- EXPRESSION SCALE CHECK ---")
print(expression.stack().describe(percentiles=[0.01, 0.25, 0.5, 0.75, 0.99]))
print("Number of negative values:", (expression < 0).sum().sum())
print("Number of zero values:", (expression == 0).sum().sum())