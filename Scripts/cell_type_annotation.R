
setwd('D:/SC-gastric/')
library(Seurat)
library(SingleR)
library(readr)
library(tidyverse)
library(celldex)


#read .rds file
gc_merged <- readRDS("gc_merged_with_clustering.rds")
gc_merged

Idents(gc_merged)  # Should show clusters like 0,1,2,...
view(gc_merged@meta.data)


Idents(gc_merged) <- "seurat_clusters"  # Set the identity class to seurat_clusters
gc_merged <- FindNeighbors(gc_merged, dims = 1:18)
gc_merged <- FindClusters(gc_merged, resolution = 0.15)  # Try 0.2 first
table(gc_merged$seurat_clusters)

# Step 2: Prepare for DE analysis with SCT
gc_merged <- PrepSCTFindMarkers(gc_merged)


markers <- FindAllMarkers(
  object = gc_merged,
  assay = "SCT",           # Use SCT assay for DE test
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)
markers
view(markers)




library(dplyr)

top20 <- markers %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC)

view(top20)

#save the file
write_csv(top20, "top20_markers_8clusters.csv")







# Extract only unique gene names
top_genes <- unique(top20$gene)
DoHeatmap(gc_merged, features = top_genes, group.by = "seurat_clusters") +
  NoLegend()


#save the plot
ggsave('Heatmap_top_20', height = 8, width = 9)









#########################################################
ref <- celldex::HumanPrimaryCellAtlasData()
view(colData(ref))

#extract test data from seurat object
test.data <- GetAssayData(gc_merged, slot = 'data')


# Run SingleR
pred <- SingleR(test = as.matrix(test.data), 
                ref = ref, 
                labels = ref$label.main)

pred

#add predicted level to seurat metadata
gc_merged$SingleR.labels <- pred$labels
#visualize annotation in UMAP or tSNE
DimPlot(gc_merged, group.by = "SingleR.labels", label = TRUE, reduction = "umap") + 
  ggtitle("Cell Type Annotation (SingleR)")

##########################################################





#inspect each object from seurat object to get a feel
view(GetAssayData(gc_merged, assay = "SCT", slot = "counts"))
view(GetAssayData(gc_merged, assay = "SCT", slot = "counts"))
Embeddings(gc_merged, reduction = "pca")[1:5, 1:5]
head(VariableFeatures(gc_merged), 10) 
view(gc_merged@meta.data)
DefaultAssay(gc_merged)


rownames(gc_merged@meta.data)[10]


#show first 5 value of cell barcode of cell_labels
head(cell_labels$cell_barcode, 5)





#attach final celltype to seurat metadata
# Extract full Seurat barcodes
seurat_barcodes <- colnames(gc_merged)

# Extract 16bp core from Seurat barcodes
seurat_core <- sapply(seurat_barcodes, function(x) {
  parts <- strsplit(x, "_")[[1]]
  sub("-\\d$", "", parts[3])
})

# Add to metadata for mapping
gc_merged$barcode_core <- seurat_core  # temporary core barcode column

# Extract 16bp core from cell_labels
cell_labels$barcode_core <- sub("-\\d$", "", cell_labels$cell_barcode)

# Now merge the two by barcode_core
library(dplyr)
metadata_annotated <- gc_merged@meta.data %>%
  left_join(cell_labels[, c("barcode_core", "final_celltype")], by = "barcode_core")

# Confirm dimensions
stopifnot(nrow(metadata_annotated) == ncol(gc_merged))

# Update Seurat metadata
gc_merged@meta.data <- metadata_annotated

view(gc_merged@meta.data)


sum(duplicated(cell_labels_unique$barcode_core))



# Show duplicated barcode_core rows
cell_labels %>%
  filter(duplicated(barcode_core) | duplicated(barcode_core, fromLast = TRUE)) %>%
  arrange(barcode_core)

view(cell_labels)
# Filter to tumor-only cells
cell_labels_tumor <- cell_labels %>%
  filter(condition == "tumor")
view(cell_labels_tumor)
#remove duplicates
cell_labels_unique <- cell_labels_tumor %>%
  distinct(barcode_core, .keep_all = TRUE)

view(cell_labels_unique)


#see how many barcode_core matches
matched_barcodes <- intersect(gc_merged@meta.data$barcode_core, cell_labels_unique$barcode_core)
length(matched_barcodes)


table(gc_merged$sample)  # assuming you renamed orig.ident to sample
table(cell_labels_unique$orig.ident)


unique(gc_merged$sample)


# Keep only annotations from samples that exist in gc_merged
valid_samples <- unique(gc_merged$sample)

cell_labels_filtered <- cell_labels_unique %>%
  filter(orig.ident %in% valid_samples)
view(cell_labels_filtered)
nrow(cell_labels_filtered)
nrow(gc_merged@meta.data)


head(colnames(gc_merged), 5)
head(cell_labels_filtered$cell_barcode, 5)





# Safe join: only matched barcodes will get annotated
metadata_annotated <- gc_merged@meta.data %>%
  left_join(cell_labels_filtered[, c("barcode_core", "final_celltype")], by = "barcode_core")

# Check if all rows are preserved
stopifnot(nrow(metadata_annotated) == ncol(gc_merged))

# Assign new metadata to Seurat object
gc_merged@meta.data <- metadata_annotated
view(gc_merged@meta.data)



# Subset Seurat object to only annotated cells
gc_annotated <- subset(gc_merged, cells = which(!is.na(gc_merged$final_celltype)))

# Now plot only annotated cells
DimPlot(gc_annotated, reduction = "umap", group.by = "final_celltype", label = TRUE, pt.size = 0.5) +
  ggtitle("UMAP: Annotated by Final Cell Type")





annotated_cells <- rownames(gc_merged@meta.data[!is.na(gc_merged$final_celltype), ])
length(annotated_cells)  # Check how many annotated cells
gc_annotated <- gc_merged[, annotated_cells]

# Clear broken graph slots (optional but good practice)
gc_merged@graphs <- list()
gc_merged@neighbors <- list()



annotated_cells <- rownames(gc_merged@meta.data[!is.na(gc_merged$final_celltype), ])
gc_annotated <- gc_merged[, annotated_cells]


head(colnames(gc_merged))
head(annotated_cells)





