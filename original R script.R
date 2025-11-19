
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(dplyr)
library(Seurat) # v5.3
library(patchwork)

# 1 Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "/Users/outobi/Documents/Learn bioinformatics /scRNA by Seurat V5/filtered_gene_bc_matrices/hg19/")

# input genes.txt.  barcodes.txt and matrix.mtx
# output pbmc.data as sparse matrix object dgCMatrix.   gene are rows and cell are columns 




# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc

## An object of class Seurat 
## 13714 features across 2700 samples within 1 assay 
## Active assay: RNA (13714 features, 0 variable features)
##  1 layer present: counts

pbmc.data[c("CD3D", "TCL1A", "MS4A1"), 1:30]


# 2 QC
# The [[ operator can add columns to object metadata directly! This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc@assays
pbmc@meta.data
#pbmc@meta.data$percent.mt = pbmc[["percent.mt"]]

  
# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# real filter
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)



# 3 normalization

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)


# 4 find variable features
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2


# 5 scale

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)


# data structure in V5
pbmc@assays$RNA@layers$counts
pbmc@assays$RNA@layers$data
pbmc@assays$RNA@layers$scale.data[c(1:5),c(1:5)]
# same as pbmc[["RNA"]]$scale.data


# 6 PCA 
# only in selected most variable features after scale
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)


#pbmc
#├── assays
#│   └── RNA
#├── meta.data
#├── reductions
#│   ├── pca
#│   ├── umap
#│   └── tsne
#└── ...

#cell.embedding.  PCA Score.  2638 cells * 50 PC.  Usigma or Z
#feature.loading. PCA loading. 2000 genes * 50 PC V 

# pc1 and pc2
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
# pc1 
VizDimLoadings(pbmc, dims = 1, reduction = "pca")

# scatter plot
DimPlot(pbmc, reduction = "pca") + NoLegend()

# pc1 component, top and bottom 250 cells along pc1 component, heatmap of scaled expression 
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)


# pc number selection
ElbowPlot(pbmc)
#We chose 10 here, but encourage users to consider the following:
  
#  Dendritic cell and NK aficionados may recognize that genes strongly associated with PCs 12 and 13 define rare immune subsets (i.e. MZB1 is a marker for plasmacytoid DCs). However, these groups are so rare, they are difficult to distinguish from background noise for a dataset of this size without prior knowledge.
#We encourage users to repeat downstream analyses with a different number of PCs (10, 15, or even 50!). As you will observe, the results often do not differ dramatically.
#We advise users to err on the higher side when choosing this parameter. For example, performing downstream analyses with only 5 PCs does significantly and adversely affect results.


# 7 cluster the cells

#•	FindNeighbors()
#•	用 PCA 前 10 个主成分做距离
#•	每个细胞找 k 近邻 → KNN graph
#•	用 Jaccard 相似度把边重新加权 → SNN graph (pbmc@graphs$RNA_snn)
#•	FindClusters()
#•	在这个 SNN 图上跑 Louvain/SLM
#•	根据 modularity 找社区 → cluster
#•	resolution 控制“粗/细”：值越大，cluster 越多
#•	结果写到 meta.data$seurat_clusters (new column) 和 Idents(pbmc)

pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)


names(pbmc@graphs)
# "RNA_nn"  "RNA_snn"
pbmc@graphs$RNA_snn   # sparse matrix of edge weights


RNA_nn <- pbmc@graphs$RNA_nn # 2638 * 2638 0 or 1 matrix (1 is in this cell k closest distance)
dim(RNA_nn)               # n_cells x n_cells
RNA_nn[1:5, 1:5]  

RNA_snn <- pbmc@graphs$RNA_snn # 2638 * 2638 0 to 1 matrix (quantify their similarity or overlap in knn)
dim(RNA_snn)
RNA_snn[1:10, 1:10]
range(RNA_snn[RNA_snn != 0])
# 0.08108108 1.00000000


# storage places
head(Idents(pbmc), 5)

pbmc$seurat_clusters


# 8 UMAP or tSNE

# input first ten PC in each cell, output 2 umap components in each cell
pbmc <- RunUMAP(pbmc, dims = 1:10)
# output store in pbmc@reductions[["umap"]]



# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters

DimPlot(pbmc, reduction = "umap", label= TRUE)


# 9 DE analysis


# find all markers of cluster 2
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2)
head(cluster2.markers, n = 5)
##             p_val avg_log2FC pct.1 pct.2    p_val_adj
## IL32 2.593535e-91  1.3221171 0.949 0.466 3.556774e-87
## LTB  7.994465e-87  1.3450377 0.981 0.644 1.096361e-82
## CD3D 3.922451e-70  1.0562099 0.922 0.433 5.379250e-66
## IL7R 1.130870e-66  1.4256944 0.748 0.327 1.550876e-62
## LDHB 4.082189e-65  0.9765875 0.953 0.614 5.598314e-61
# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3))
head(cluster5.markers, n = 5)
##                       p_val avg_log2FC pct.1 pct.2     p_val_adj
## FCGR3A        2.150929e-209   6.832372 0.975 0.039 2.949784e-205
## IFITM3        6.103366e-199   6.181000 0.975 0.048 8.370156e-195
## CFD           8.891428e-198   6.052575 0.938 0.037 1.219370e-193
## CD68          2.374425e-194   5.493138 0.926 0.035 3.256286e-190
## RP11-290F20.3 9.308287e-191   6.335402 0.840 0.016 1.276538e-186
# find markers for every cluster compared to all remaining cells, report only the positive
# ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)
pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)
## # A tibble: 7,046 × 7
## # Groups:   cluster [9]
##        p_val avg_log2FC pct.1 pct.2 p_val_adj cluster gene     
##        <dbl>      <dbl> <dbl> <dbl>     <dbl> <fct>   <chr>    
##  1 1.74e-109       1.19 0.897 0.593 2.39e-105 0       LDHB     
##  2 1.17e- 83       2.37 0.435 0.108 1.60e- 79 0       CCR7     
##  3 8.94e- 79       1.09 0.838 0.403 1.23e- 74 0       CD3D     
##  4 3.05e- 53       1.02 0.722 0.399 4.19e- 49 0       CD3E     
##  5 3.28e- 49       2.10 0.333 0.103 4.50e- 45 0       LEF1     
##  6 6.66e- 49       1.25 0.623 0.358 9.13e- 45 0       NOSIP    
##  7 9.31e- 44       2.02 0.328 0.11  1.28e- 39 0       PRKCQ-AS1
##  8 4.69e- 43       1.53 0.435 0.184 6.43e- 39 0       PIK3IP1  
##  9 1.47e- 39       2.70 0.195 0.04  2.01e- 35 0       FHIT     
## 10 2.44e- 33       1.94 0.262 0.087 3.34e- 29 0       MAL      
## # ℹ 7,036 more rows


#Seurat has several tests for differential expression which can be set with the test.use parameter (see our DE vignette for details). For example, the ROC test returns the ‘classification power’ for any individual marker (ranging from 0 - random, to 1 - perfect).

cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)


# classical visualization
VlnPlot(pbmc, features = c("MS4A1", "CD79A"))

# you can plot raw counts as well
VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)


FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
                               "CD8A"))

#1.	group_by(cluster)
#→ “Split” the table by cluster.
#Now all subsequent operations happen within each cluster separately.
#2.	filter(avg_log2FC > 1)
#→ Keep only genes that are strongly upregulated in that cluster
#(log2 fold change > 1 means at least ~2× higher expression).
#3.	slice_head(n = 10)
#→ For each cluster, keep the first 10 rows (after filtering).
#Since FindAllMarkers usually sorts by p_val or p_val_adj, this means:
#  “Take the top 10 strongest / most significant marker genes per cluster.”
#4.	ungroup()
#→ Remove the grouping so top10 is now just a normal tibble.
#5.	-> top10
#→ Save the result as a new object top10.
pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10

#Draw a heatmap of the top 10 upregulated marker genes for each cluster across all cells.
DoHeatmap(pbmc, features = top10$gene) + NoLegend()




# Single gene
RidgePlot(pbmc, features = "MS4A1")

# Multiple genes, faceted
RidgePlot(pbmc, features = c("MS4A1", "CD79A", "CD3D"), ncol = 1)


# feature correlation among all cells
FeatureScatter(pbmc, feature1 = "MS4A1", feature2 = "CD79A")

FeatureScatter(pbmc, feature1 = "CD3D", feature2 = "CCR7")

# cell correlation among variable features
CellScatter(
  pbmc,
  cell1 = colnames(pbmc)[1],
  cell2 = colnames(pbmc)[2],
  features = VariableFeatures(pbmc)
)


# classical bubble plot of marker gene
markers.to.plot <- c("MS4A1", "CD79A",    # B cells
                     "CD3D", "IL7R",      # CD4 T
                     "LTB", "NKG7")       # others

DotPlot(pbmc, features = markers.to.plot) + RotatedAxis()



# 10 cell type labeling
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()


saveRDS(pbmc, file = "../pbmc_tutorial.rds") # save r data serialization, one object only
