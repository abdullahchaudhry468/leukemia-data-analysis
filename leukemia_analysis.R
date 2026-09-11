library(golubEsets)
library(Biobase)
library(vsn)
library(limma)
library(ggplot2)
library(pheatmap)
library(AnnotationDbi)
library(hu6800.db)

data(Golub_Merge)

project_dir <- getwd()
data_dir <- file.path(project_dir, "data")
figures_dir <- file.path(project_dir, "figures")

dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

expression_data <- exprs(Golub_Merge)
sample_data <- pData(Golub_Merge)

stopifnot(
  identical(colnames(expression_data), rownames(sample_data)),
  sum(is.na(expression_data)) == 0
)

normalized_data <- justvsn(Golub_Merge)
expr_norm <- exprs(normalized_data)

write.csv(
  expr_norm,
  file.path(data_dir, "expression_normalized.csv")
)

write.csv(
  sample_data,
  file.path(data_dir, "sample_metadata.csv")
)

png(
  file.path(figures_dir, "vsn_mean_sd_plot.png"),
  width = 1200,
  height = 900,
  res = 150
)

meanSdPlot(normalized_data)
dev.off()

group <- factor(pData(normalized_data)$ALL.AML)

stopifnot(
  identical(colnames(expr_norm), rownames(pData(normalized_data)))
)

design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

fit <- lmFit(expr_norm, design)

contrast_matrix <- makeContrasts(
  AMLvsALL = AML - ALL,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

results <- topTable(
  fit2,
  coef = "AMLvsALL",
  number = Inf,
  adjust.method = "BH"
)

gene_symbols <- AnnotationDbi::mapIds(
  hu6800.db,
  keys = rownames(results),
  column = "SYMBOL",
  keytype = "PROBEID",
  multiVals = "first"
)

results$GeneSymbol <- gene_symbols[rownames(results)]

significant <- results[
  results$adj.P.Val < 0.05 &
    abs(results$logFC) > 1,
]

write.csv(
  results,
  file.path(data_dir, "differential_expression_annotated.csv")
)

write.csv(
  significant,
  file.path(data_dir, "differential_expression_significant_annotated.csv")
)

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

ggsave(
  file.path(figures_dir, "volcano_plot.png"),
  plot = volcano_plot,
  width = 8,
  height = 6,
  dpi = 300
)

top_probes <- rownames(results)[
  order(results$adj.P.Val)
][1:50]

heatmap_matrix <- expr_norm[top_probes, ]

annotation_col <- data.frame(
  Leukemia = group
)

rownames(annotation_col) <- colnames(heatmap_matrix)

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
  filename = file.path(figures_dir, "heatmap_top50.png"),
  width = 9,
  height = 7
)

annotated_significant <- significant[
  !is.na(significant$GeneSymbol),
]

annotated_significant <- annotated_significant[
  order(annotated_significant$adj.P.Val),
]

top20_genes <- annotated_significant[1:20, ]

write.csv(
  top20_genes,
  file.path(data_dir, "top20_annotated_genes.csv")
)
