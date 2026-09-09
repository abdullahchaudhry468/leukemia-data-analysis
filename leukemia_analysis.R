# ============================================================
# Leukemia Gene Expression Analysis
# Golub dataset: ALL vs AML
# ============================================================

# Load the packages used in this analysis
library(golubEsets)
library(Biobase)
library(vsn)
library(limma)
library(ggplot2)
library(pheatmap)
library(AnnotationDbi)
library(hu6800.db)

# Load the published leukemia dataset
data(Golub_Merge)

# Set up the project folders
project_dir <- "~/Documents/GitHub/leukemia-bioinformatics"
data_dir <- file.path(project_dir, "data")
figures_dir <- file.path(project_dir, "figures")

dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# Check the dataset
print(Golub_Merge)
# ============================================================
# 1. Extract and validate the original data
# ============================================================

# Extract expression measurements and sample metadata
expression_data <- exprs(Golub_Merge)
sample_data <- pData(Golub_Merge)

# Check dimensions and leukemia groups
print(dim(expression_data))
print(dim(sample_data))
print(table(sample_data$ALL.AML))

# Confirm that sample identifiers are aligned
stopifnot(
  identical(colnames(expression_data), rownames(sample_data))
)

# Check for missing expression measurements
print(sum(is.na(expression_data)))

# Save the original data
write.csv(
  expression_data,
  file.path(data_dir, "expression.csv")
)

write.csv(
  sample_data,
  file.path(data_dir, "sample_metadata.csv")
)
# ============================================================
# 2. Normalize the microarray data and perform QC
# ============================================================

# Apply variance-stabilizing normalization
normalized_data <- justvsn(Golub_Merge)

# Extract normalized expression measurements
expr_norm <- exprs(normalized_data)

# Save normalized data
write.csv(
  expr_norm,
  file.path(data_dir, "expression_normalized.csv")
)

saveRDS(
  normalized_data,
  file.path(data_dir, "normalized_golub.rds")
)

# Inspect normalized expression values
print(summary(as.vector(expr_norm)))

# Save the mean-SD diagnostic plot
png(
  filename = file.path(figures_dir, "vsn_mean_sd_plot.png"),
  width = 1200,
  height = 900,
  res = 150
)

meanSdPlot(normalized_data)

dev.off()

# ============================================================
# 3. Differential expression: AML vs ALL
# ============================================================

# Create leukemia group labels
group <- factor(pData(normalized_data)$ALL.AML)

# Confirm sample alignment
stopifnot(
  identical(colnames(expr_norm), rownames(pData(normalized_data)))
)

print(table(group))

# Create the design matrix
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

# Fit the linear model
fit <- lmFit(expr_norm, design)

# Define the contrast: AML minus ALL
contrast_matrix <- makeContrasts(
  AMLvsALL = AML - ALL,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# Extract all differential-expression results
results <- topTable(
  fit2,
  coef = "AMLvsALL",
  number = Inf,
  adjust.method = "BH"
)

# Select significant probe sets
significant <- results[
  results$adj.P.Val < 0.05 &
    abs(results$logFC) > 1,
]

# Display the number of significant probes
print(nrow(significant))
print(sum(significant$logFC > 0))
print(sum(significant$logFC < 0))

# Save the results
write.csv(
  results,
  file.path(data_dir, "differential_expression_all.csv")
)

write.csv(
  significant,
  file.path(data_dir, "differential_expression_significant.csv")
)
# ============================================================
# 4. Volcano plot
# ============================================================

# Classify probe sets using the significance thresholds
volcano_data <- results
volcano_data$Significance <- "Not significant"

volcano_data$Significance[
  volcano_data$adj.P.Val < 0.05 &
    volcano_data$logFC > 1
] <- "Higher in AML"

volcano_data$Significance[
  volcano_data$adj.P.Val < 0.05 &
    volcano_data$logFC < -1
] <- "Higher in ALL"

# Create the volcano plot
volcano_plot <- ggplot(
  volcano_data,
  aes(x = logFC, y = -log10(adj.P.Val), color = Significance)
) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c(
    "Higher in AML" = "#C75B5B",
    "Higher in ALL" = "#4C78A8",
    "Not significant" = "grey70"
  )) +
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) +
  labs(
    title = "Differential Expression: AML vs ALL",
    x = "Log2 fold change (AML vs ALL)",
    y = "-Log10 adjusted p-value",
    color = "Classification"
  ) +
  theme_minimal()

# Save the plot
ggsave(
  filename = file.path(figures_dir, "volcano_plot.png"),
  plot = volcano_plot,
  width = 8,
  height = 6,
  dpi = 300
)
# ============================================================
# 5. Heatmap of top differentially expressed probes
# ============================================================

# Select the 50 probes with the lowest adjusted p-values
top_probes <- rownames(results)[
  order(results$adj.P.Val)
][1:50]

# Create expression matrix for those probes
heatmap_matrix <- expr_norm[top_probes, ]

# Create sample annotation
annotation_col <- data.frame(
  Leukemia = group
)

rownames(annotation_col) <- colnames(heatmap_matrix)

# Generate and save the heatmap
pheatmap(
  heatmap_matrix,
  scale = "row",
  annotation_col = annotation_col,
  show_colnames = FALSE,
  show_rownames = FALSE,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  main = "Top 50 Differentially Expressed Probes",
  filename = file.path(
    figures_dir,
    "heatmap_top50.png"
  ),
  width = 9,
  height = 7
)
# ============================================================
# 6. Annotate significant probes with gene symbols
# ============================================================

# Map probe IDs to gene symbols
gene_symbols <- AnnotationDbi::mapIds(
  hu6800.db,
  keys = rownames(results),
  column = "SYMBOL",
  keytype = "PROBEID",
  multiVals = "first"
)

# Add gene symbols to the full results table
results$GeneSymbol <- gene_symbols[rownames(results)]

# Recreate the significant-results table with annotation
significant <- results[
  results$adj.P.Val < 0.05 &
    abs(results$logFC) > 1,
]

# Save annotated results
write.csv(
  results,
  file.path(data_dir, "differential_expression_annotated.csv")
)

write.csv(
  significant,
  file.path(data_dir, "differential_expression_significant_annotated.csv")
)

# Keep only significant probes with mapped gene symbols
annotated_significant <- significant[
  !is.na(significant$GeneSymbol),
]

# Order by adjusted p-value
annotated_significant <- annotated_significant[
  order(annotated_significant$adj.P.Val),
]

# Select the top 20 annotated probes
top20_genes <- annotated_significant[1:20, ]

print(top20_genes)

# Save the top 20 table
write.csv(
  top20_genes,
  file.path(data_dir, "top20_annotated_genes.csv")
)
