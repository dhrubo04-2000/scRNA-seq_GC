# Load necessary libraries
library(Seurat)
library(SeuratObject)
library(ggplot2)
library(sctransform)
library(glmGamPoi)
setwd('D:/SC-gastric/gastric_scRNAseq_raw/')
# Define all sample paths (correct duplicated key "5931_n1")
sample_paths <- list(
  "5846_t1" = "D:/SC-gastric/gastric_scRNAseq_raw/5846_t1",
  "5866_t1" = "D:/SC-gastric/gastric_scRNAseq_raw/5866_t1",
  "5931_t1" = "D:/SC-gastric/gastric_scRNAseq_raw/5931_t1",
  "6207_t1" = "D:/SC-gastric/gastric_scRNAseq_raw/6207_t1",
  "6592_t1" = "D:/SC-gastric/gastric_scRNAseq_raw/6592_t1",
  "6709_t1" = "D:/SC-gastric/gastric_scRNAseq_raw/6709_t1",
  "6649_t1" = "D:/SC-gastric/gastric_scRNAseq_raw/6649_t1")

# Create Seurat objects and store in list
seurat_list <- list()

for (sample_name in names(sample_paths)) {
  path <- sample_paths[[sample_name]]
  data <- Read10X(data.dir = path)
  
  seu <- CreateSeuratObject(counts = data, project = sample_name,
                            min.cells = 3, min.features = 200)
  
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
  
  seu <- subset(seu, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20)
  
  seu <- RenameCells(seu, add.cell.id = sample_name)
  seu$sample <- sample_name
  seu$group <- ifelse(grepl("_n", sample_name), "Normal", "Tumor")  # Add group info
  
  seurat_list[[sample_name]] <- seu
}

# Merge all objects
gc_merged <- Reduce(function(x, y) merge(x, y, project = "GC_all"), seurat_list)

#save gc_merged as .rds
saveRDS(gc_merged, file = "gc_merged.rds")

# Check number of cells per sample and per group
table(gc_merged$sample)
table(gc_merged$group)


VlnPlot(gc_merged, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        group.by = "sample", 
        pt.size = 0.1, 
        ncol = 3,
        layer = "counts")  # raw counts

#save with ggsave
ggsave("violin_plot_features.png", width = 12, height = 6)


# Scatter plot: Total counts vs percent mitochondrial reads
plot1 <- FeatureScatter(gc_merged, feature1 = "nCount_RNA", feature2 = "percent.mt")

# Scatter plot: Total counts vs number of detected features (genes)
plot2 <- FeatureScatter(gc_merged, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# Combine plots side by side
plot1 + plot2






library(ggplot2)

# Extract metadata as data frame
meta_df <- gc_merged@meta.data

# Histogram: Total counts
ggplot(meta_df, aes(x = nCount_RNA)) +
  geom_histogram(bins = 80, fill = "skyblue", color = "black") +
  ggtitle("Histogram of Total Counts per Cell") +
  xlab("Total Counts (nCount_RNA)") + ylab("Frequency")

# Histogram: Number of genes
ggplot(meta_df, aes(x = nFeature_RNA)) +
  geom_histogram(bins = 80, fill = "lightgreen", color = "black") +
  ggtitle("Histogram of Number of Genes per Cell") +
  xlab("Number of Genes Detected (nFeature_RNA)") + ylab("Frequency")



library(ggplot2)

# Extract total counts per cell (barcode)
counts <- gc_merged$nCount_RNA

# Sort counts decreasingly and get ranks
sorted_counts <- sort(counts, decreasing = TRUE)
ranks <- seq_along(sorted_counts)

# Create data frame for plotting
df <- data.frame(rank = ranks, count = sorted_counts)

# Plot barcode rank plot (log-log scale)
ggplot(df, aes(x = rank, y = count)) +
  geom_line(color = "blue") +
  scale_x_log10() +
  scale_y_log10() +
  coord_cartesian(ylim = c(1e2, 1e5)) +  # zoom out y-axis
  ggtitle("Barcode Rank Plot") +
  xlab("Barcode Rank (log scale)") +
  ylab("Total Counts per Cell (log scale)") +
  theme_minimal()













# Normalize with SCTransform
gc_merged <- SCTransform(gc_merged, vars.to.regress = "percent.mt", 
                         verbose = FALSE, method = "glmGamPoi")



str(gc_merged)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(gc_merged), 10)
top10














# PCA
gc_merged <- RunPCA(gc_merged, verbose = FALSE)
ElbowPlot(gc_merged)  # Use this to pick dimensions




# Clustering & UMAP using top PCs (adjust if needed)
gc_merged <- FindNeighbors(gc_merged, dims = 1:17)
gc_merged <- FindClusters(gc_merged, resolution = 1.0)

#visualize dimplot for best resolution
DimPlot(gc_merged, group.by = "RNA_snn_res.1.0", label = TRUE)
DimPlot(gc_merged, label = TRUE)


Idents(gc_merged) <- "RNA_snn_res.1.0"
Idents(gc_merged) 
#run tSNE
gc_merged <- RunTSNE(gc_merged, dims = 1:17)
#run umap
gc_merged <- RunUMAP(gc_merged, dims = 1:17)

#save gc_mergeed as .rds
saveRDS(gc_merged, file = "gc_merged_with_clustering.rds")


# Plot tSNE with clusters
DimPlot(gc_merged, reduction = "tsne", group.by = "seurat_clusters", label = TRUE, pt.size = 0.5) +
  ggtitle("t-SNE Plot of Clusters")
# Plot UMAP with clusters
DimPlot(gc_merged, reduction = "umap", group.by = "seurat_clusters", label = TRUE, pt.size = 0.5) +
  ggtitle("UMAP Plot of Clusters")

#save plot using ggsave
ggsave("tsne_plot_clusters.png", width = 8, height = 6)




                              


















