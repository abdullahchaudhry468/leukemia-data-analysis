# leukemia-bioinformatics
Reproducible analysis of leukemia gene-expression data using R, Python, statistical learning and bioinformatics workflows.

## Project Overview

This project explores gene-expression differences between acute lymphoblastic leukemia (ALL) and acute myeloid leukemia (AML) using the publicly available Golub leukemia microarray dataset. The goal is to develop a reproducible bioinformatics workflow and investigate how expression patterns differ between the two leukemia types.

The analysis combines R and Python to perform data validation, variance-stabilizing normalization, principal component analysis (PCA), differential-expression testing, and visualization. Gene symbols are mapped to significant probe sets to support biological interpretation.

This is an educational and exploratory analysis of a historical dataset. It is not intended for clinical diagnosis or to establish new biomarkers.

## Dataset

This project uses the Golub leukemia gene-expression dataset, originally published by Golub et al. (1999) and distributed through the Bioconductor `golubEsets` package.

The dataset contains 7,129 microarray probe sets measured across 72 leukemia samples: 47 acute lymphoblastic leukemia (ALL) samples and 25 acute myeloid leukemia (AML) samples. The data are historical and have undergone preprocessing prior to distribution.

The analysis uses `Golub_Merge`, which combines the original training and test samples. The dataset is used for exploratory analysis and differential-expression testing, not for clinical prediction.

**Original publication:** Golub TR et al. (1999). *Molecular Classification of Cancer: Class Discovery and Class Prediction by Gene Expression Monitoring*. Science, 286, 531–537.

**Data source:** https://bioconductor.org/packages/golubEsets/


## Methods

### Data validation and normalization

The expression matrix and sample metadata were extracted from `Golub_Merge`. Sample identifiers, group counts, and missing values were checked before analysis. Variance-stabilizing normalization was performed using the Bioconductor `vsn` package, and a mean–SD diagnostic plot was generated to assess the normalized data.

### Principal component analysis

Python was used to perform principal component analysis (PCA) on the expression data. Features were standardized using `StandardScaler`, and the first two principal components were visualized to explore variation and sample grouping between ALL and AML.

### Differential expression

Differential-expression analysis was performed in R using `limma`. A linear model compared AML against ALL, followed by empirical Bayes moderation. P-values were adjusted using the Benjamini–Hochberg method.

Probe sets were considered significant when they met both criteria:

* Adjusted p-value < 0.05
* Absolute log2 fold change > 1

Positive log2 fold changes indicate higher expression in AML, while negative values indicate higher expression in ALL.

### Visualization and annotation

A volcano plot was created to display differential-expression results. A heatmap of the 50 probe sets with the lowest adjusted p-values was generated using row-scaled expression values and hierarchical clustering.

Significant probe sets were mapped to gene symbols using `hu6800.db` and `AnnotationDbi`. The top 20 annotated probe sets were exported for further biological interpretation. Probe sets are not necessarily unique genes, and some historical probe annotations may be ambiguous or unmapped.


## Results

### Differential-expression summary

Using an adjusted p-value threshold of 0.05 and an absolute log2 fold-change threshold of 1, the analysis identified **328 significant probe sets**:

* 170 probe sets with higher expression in AML.
* 158 probe sets with higher expression in ALL.

Of the 328 significant probe sets, 314 were mapped to gene symbols using the historical microarray annotation package, while 14 remained unmapped. These counts represent probe sets rather than necessarily unique genes.

### Principal component analysis

The PCA visualization showed partial separation between ALL and AML samples, with some overlap and outliers. The first two principal components explained approximately 14.9% and 9.4% of the variance, respectively. This suggests that leukemia subtype contributes to expression variation, although the first two components do not completely separate the groups.

![PCA of leukemia samples](figures/pca_leukemia.png)

### Volcano plot

The volcano plot displays the magnitude and statistical significance of expression differences between AML and ALL. Positive log2 fold changes indicate higher expression in AML, while negative values indicate higher expression in ALL.

![Differential-expression volcano plot](figures/volcano_plot.png)

### Heatmap

The heatmap of the top 50 differentially expressed probe sets shows broad clustering by leukemia subtype, with some sample mixing. Expression values were scaled by row, so the colors represent relative expression within each probe set rather than absolute expression levels.

![Top 50 differentially expressed probes](figures/heatmap_top50.png)

### Selected annotated findings

Among the highly ranked annotated probe sets were CD33, MPO, CFD, CST3, DNTT, and TCF3. These results provide candidates for further biological interpretation, but their presence in this exploratory analysis does not establish them as novel or clinically validated biomarkers.


## Reproducibility

### Requirements

The analysis uses R and Python. The main R dependencies are `golubEsets`, `Biobase`, `vsn`, `limma`, `ggplot2`, `pheatmap`, `AnnotationDbi`, and `hu6800.db`. The Python workflow uses pandas, NumPy, scikit-learn, and Matplotlib.

### Running the analysis

1. Clone or download this repository.
2. Install the required R and Python packages.
3. Open `leukemia_analysis.R` in RStudio and update `project_dir` to match the location of the repository on your computer.
4. Run the R script from beginning to end. It loads the dataset, exports the data, performs normalization and differential-expression analysis, and generates the R figures and result tables.
5. Run `leukemia_analysis.py` to reproduce the Python PCA analysis.

The R script saves output files to the `data/` and `figures/` folders. The Python script also saves its PCA figure to `figures/`.

**Note:** The Python PCA workflow and the R differential-expression workflow should be interpreted as separate exploratory analyses unless the same normalized input is explicitly used in both. The R workflow uses VSN-normalized expression values.

## Limitations

This project uses a historical microarray dataset with 72 samples, and the findings may not generalize to modern or independent patient cohorts. The data were previously processed before distribution, and additional normalization was applied for this exploratory workflow.

Differential-expression results are reported at the probe-set level. Multiple probes may correspond to the same gene, and historical annotations can be incomplete or ambiguous. The heatmap uses the same group labels that were used to select the top probes, so its apparent separation is descriptive rather than independent validation.

No independent clinical validation, diagnostic model, or prospective testing was performed. The results should therefore be considered exploratory and educational, not evidence of new biomarkers or clinical diagnostic performance.


## References

1. Golub TR, Slonim DK, Tamayo P, et al. (1999). Molecular Classification of Cancer: Class Discovery and Class Prediction by Gene Expression Monitoring. *Science*, 286, 531–537. https://doi.org/10.1126/science.286.5439.531

2. Bioconductor. `golubEsets`: Golub leukemia data. https://bioconductor.org/packages/golubEsets/

3. Ritchie ME, Phipson B, Wu D, et al. (2015). limma powers differential expression analyses for RNA-sequencing and microarray studies. *Nucleic Acids Research*, 43(7), e47. https://doi.org/10.1093/nar/gkv007

4. Huber W, von Heydebreck A, Sültmann H, Poustka A, Vingron M. (2002). Variance stabilization applied to microarray data calibration and to the quantification of differential expression. *Bioinformatics*, 18(Suppl 1), S96–S104. https://doi.org/10.1093/bioinformatics/18.suppl_1.S96

## Project Status

The exploratory analysis and primary visualizations are complete. The repository is being prepared for public release, including final checks of reproducibility, documentation, and file organization.


## PCA Analysis

Principal component analysis (PCA) was used to explore overall patterns in gene-expression profiles across the 72 leukemia samples. The first two principal components explained 10.5% and 5.0% of the total variance, respectively. The PCA showed clear separation between ALL and AML samples, indicating distinct global gene-expression patterns between the two leukemia subtypes.

![PCA of leukemia gene expression](figures/pca_leukemia.png)

## Analysis Workflow

1. Loaded the Golub leukemia gene-expression dataset and corresponding sample metadata.
2. Performed variance-stabilizing normalization and quality checks in R.
3. Compared gene expression between ALL and AML samples using differential-expression analysis.
4. Applied multiple-testing correction using the Benjamini-Hochberg method.
5. Identified significantly differentially expressed probes using adjusted p-value < 0.05 and |log2 fold change| > 1.
6. Annotated significant probes with gene symbols using the `hu6800.db` annotation package.
7. Visualized differential expression using a volcano plot and hierarchical clustering heatmap.
8. Used PCA to examine global expression patterns and assess separation between leukemia subtypes.
9. Exported analysis results and figures for reproducibility.

## Key Results

The analysis identified 328 significant probes meeting the selected statistical and fold-change criteria. Of these, 314 were successfully mapped to gene symbols. PCA showed separation between ALL and AML samples, with PC1 explaining 10.5% of the variance and PC2 explaining 5.0%.

## Reproducibility

The analysis code, generated figures, and exported result tables are included in this repository. The workflow is organized so that the analysis can be reproduced from the provided data and scripts.
