# Single-Cell RNA-seq Analysis with Seurat V5

A comprehensive single-cell RNA sequencing analysis workflow using the PBMC3k dataset and Seurat V5. This project demonstrates best practices for quality control, normalization, clustering, visualization, and cell type annotation of scRNA-seq data.

## 📊 Project Overview

This analysis processes 2,700 Peripheral Blood Mononuclear Cells (PBMCs) sequenced with the 10X Genomics platform. The workflow identifies 9 distinct cell populations and annotates them based on canonical marker genes.

### Identified Cell Types
- Naive CD4+ T cells
- Memory CD4+ T cells  
- CD8+ T cells
- B cells
- CD14+ Monocytes
- FCGR3A+ Monocytes
- NK cells
- Dendritic Cells
- Platelets

## 🔬 Analysis Pipeline

The complete workflow includes:

1. **Data Loading** - Import 10X Genomics data (matrix.mtx, genes.tsv, barcodes.tsv)
2. **Quality Control** - Filter cells based on gene count and mitochondrial content
3. **Normalization** - Log-normalization with scale factor of 10,000
4. **Feature Selection** - Identify 2,000 highly variable genes
5. **Scaling** - Center and scale gene expression
6. **Dimensionality Reduction** - PCA analysis
7. **Clustering** - Graph-based clustering using Louvain algorithm
8. **Visualization** - UMAP projection for 2D visualization
9. **Differential Expression** - Identify marker genes for each cluster
10. **Cell Type Annotation** - Assign biological labels based on markers

## 📋 Requirements

### Software
- R (≥ 4.0.0)
- RStudio (recommended)

### R Packages
```r
# Required packages
install.packages("dplyr")
install.packages("patchwork")

# Seurat V5
install.packages("Seurat")  # v5.3 or later
```

## 📁 Project Structure

```
filtered_gene_bc_matrices/
├── scRNA_analysis.Rmd          # Main R Markdown analysis
├── README.md                    # This file
├── hg19/                        # Data directory
│   ├── barcodes.tsv            # Cell barcodes
│   ├── genes.tsv               # Gene names
│   └── matrix.mtx              # Expression matrix
└── pbmc_tutorial.rds           # Saved Seurat object (generated)
```

## 🚀 Usage

### Running the Analysis

1. **Clone this repository:**
   ```bash
   git clone https://github.com/YOUR-USERNAME/pbmc_Seurat_scRNAseq.git
   cd pbmc_Seurat_scRNAseq
   ```

2. **Download the PBMC dataset:**
   ```r
   # The data should be placed in the hg19/ directory
   # Expected files: barcodes.tsv, genes.tsv, matrix.mtx
   ```

3. **Open and run the analysis:**
   - Open `scRNA_analysis.Rmd` in RStudio
   - Click "Knit" to generate an HTML report
   - Or run chunks interactively

4. **View results:**
   - The knitted HTML file will contain all plots and results
   - The final Seurat object is saved as `pbmc_tutorial.rds`

### Key Code Snippets

**Load data:**
```r
pbmc.data <- Read10X(data.dir = "hg19/")
pbmc <- CreateSeuratObject(counts = pbmc.data, 
                          project = "pbmc3k", 
                          min.cells = 3, 
                          min.features = 200)
```

**Quality control:**
```r
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & 
                              nFeature_RNA < 2500 & 
                              percent.mt < 5)
```

**Standard workflow:**
```r
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, nfeatures = 2000)
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc)
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:10)
```

## 📊 Key Results

### Data Quality
- **Input:** 13,714 genes × 2,700 cells
- **After QC:** 13,714 genes × 2,638 cells (62 low-quality cells removed)
- **Filtering criteria:** 200-2,500 genes per cell, <5% mitochondrial content

### Clustering
- **Method:** Louvain algorithm on SNN graph
- **Resolution:** 0.5
- **Result:** 9 distinct clusters identified

### Variable Features
- **Top 10 most variable genes:** Identified using variance stabilizing transformation
- **Used for PCA:** 2,000 highly variable genes

### Dimensionality
- **PCA:** First 10 principal components used for downstream analysis
- **UMAP:** 2D projection for visualization

## 📈 Visualizations

The analysis generates multiple publication-quality plots:

- **QC Metrics:** Violin and scatter plots for quality assessment
- **Variable Features:** Variance plots showing most informative genes
- **PCA:** Loadings, scatter plots, and heatmaps
- **Elbow Plot:** PC selection guide
- **UMAP:** Cluster visualization with cell type labels
- **Marker Expression:** Violin plots, feature plots, ridge plots
- **Heatmaps:** Top marker genes per cluster
- **Dot Plots:** Canonical marker gene expression

## 🔍 Differential Expression

Marker genes are identified using:
- Wilcoxon rank-sum test (default)
- ROC analysis for classification power
- Filters: log2FC > 1, adjusted p-value < 0.05

Top markers are visualized in heatmaps and dot plots.

## 📚 Data Source

**PBMC3k Dataset:**
- 2,700 Peripheral Blood Mononuclear Cells
- Sequenced on Illumina NextSeq 500
- 10X Genomics Chromium platform
- Human reference genome: hg19

## 🛠️ Technical Details

### Seurat Object Structure
The analysis progressively builds a comprehensive Seurat object containing:
- Raw counts, normalized data, and scaled data
- Cell-level metadata (QC metrics, cluster assignments)
- PCA and UMAP dimensionality reductions
- Graph structures (KNN and SNN)
- Marker gene lists

### Computational Considerations
- Uses sparse matrix storage for memory efficiency
- PCA computed on 2,000 variable features (not all genes)
- Scaling converts to dense matrix (may require significant memory)

## 📝 Notes

- All comments and explanations are embedded in the R Markdown file
- Code chunks include detailed input/output specifications
- Data structure and dimensions are documented at each step
- Seurat object structure is outlined after major transformations

## 🤝 Contributing

Feel free to open issues or submit pull requests for improvements.

## 📄 License

This project is open source and available under the MIT License.

## 📧 Contact

For questions or feedback, please open an issue on GitHub.

## 🙏 Acknowledgments

- **Seurat:** Stuart et al., Cell 2019
- **PBMC Dataset:** 10X Genomics
- **Tutorial inspiration:** Seurat vignettes and documentation

## 📖 References

1. Stuart T, Butler A, et al. Comprehensive Integration of Single-Cell Data. *Cell* 177, 1888-1902 (2019).
2. Hao Y, Hao S, et al. Integrated analysis of multimodal single-cell data. *Cell* 184, 3573-3587 (2021).
3. 10X Genomics PBMC Dataset: https://support.10xgenomics.com/single-cell-gene-expression/datasets

---

⭐ If you find this analysis helpful, please consider giving it a star!
