setwd('D:/SC-gastric')
library(Seurat)  
library(tidyverse)

gc_annot <- readRDS('RDS/gc_merged_annotated_final.rds')

#make umap plot of cell types
DimPlot(gc_annot, group.by = "Broad_Annotation", label = TRUE) + NoLegend()
#save the plot
ggsave('plots/umap_cell_types.png', width = 8, height = 6)


Idents(gc_annot) <- "Broad_Annotation"
gc_annot <- PrepSCTFindMarkers(gc_annot)

#find upregulated DE markers
all_markers <- FindAllMarkers(
  object = gc_annot,
  only.pos = TRUE,            # only return upregulated markers
  min.pct = 0.15,             # only test genes detected in at least 25% of cells
  logfc.threshold = 0.25      # only genes with log2 fold change > 0.25
)


view(all_markers)

#read OXPHOS_MITO_genes genes
OXPHOS_MITO_genes <- read_csv("D:/SC-gastric/OXPHOS_MITO_genes.csv")
view(OXPHOS_MITO_genes)


#intersect matching genes
matching_genes <- intersect(all_markers$gene, OXPHOS_MITO_genes$Genes)
view(matching_genes)

#rename matching genes column
matching_genes <- as.data.frame(matching_genes)
colnames(matching_genes) <- "Genes"


matching_markers <- all_markers[all_markers$gene %in% matching_genes$Genes, ]
View(matching_markers)

#save the file
write.csv(matching_markers, "matching_markers_oxphos_mito.csv", row.names = T)


###################################################################################
#read mitochondrial dynamics genes
library(readr)
Mitochondrial_Dynamics_Genes <- read_csv("Mitochondrial_Dynamics_Genes.csv")




#intersect matching genes
matching_genes <- intersect(all_markers$gene, Mitochondrial_Dynamics_Genes$Genes)


#rename matching genes column
matching_genes <- as.data.frame(matching_genes)
colnames(matching_genes) <- "Genes"


matching_markers <- all_markers[all_markers$gene %in% matching_genes$Genes, ]
View(matching_markers)








#intersect matching genes
matching_genes <- intersect(all_markers$gene, matching_genes$Genes)
view(matching_genes)

#rename matching genes column
matching_genes <- as.data.frame(matching_genes)
colnames(matching_genes) <- "Genes"


matching_markers <- all_markers[all_markers$gene %in% matching_genes$Genes, ]
View(matching_markers)


#correlation analysis of 14 matching_markers
# Extract SCT assay data (normalized expression matrix)
expr_matrix <- GetAssayData(gc_annot, assay = "SCT", slot = "data")

# Subset to only the 14 matching mitochondrial dynamics genes
matching_gene_names <- matching_genes$Genes
expr_subset <- expr_matrix[rownames(expr_matrix) %in% matching_gene_names, ]
expr_subset_t <- t(as.matrix(expr_subset))



library(corrplot)
# Pearson correlation between genes
gene_cor <- cor(expr_subset_t, method = "pearson")

corrplot(gene_cor, method = "color", type = "upper",
         tl.cex = 0.8, number.cex = 0.7, addCoef.col = "black",
         title = "Gene Correlation Heatmap", mar = c(0,0,1,0))



# Save the correlation plot
ggsave('plots/gene_correlation_heatmap_oxphos.png', width = 8, height = 6)
#save matching matching_markers
write.csv(matching_markers, "matching_markers_mitochondrial_dynamics.csv", row.names = FALSE)






#################################################################################################################
#AUC Cell calculation for Mitochondrial Oxpohos genes
library(AUCell)
# Prepare gene set list
geneSetList <- list(OXPHOS = matching_genes$Genes[1:20])
# Extract expression matrix from SCT assay (normalized data)
expr_matrix <- as.matrix(GetAssayData(gc_annot, assay = "SCT", slot = "data"))

cells_AUC <- AUCell::AUCell_buildRankings(expr_matrix, nCores=1, plotStats=FALSE)

cells_AUC_cal <- AUCell_calcAUC(geneSetList, cells_AUC, aucMaxRank = 350)


# Add AUC scores to Seurat object metadata
library(SummarizedExperiment)
#Extract the OXPHOS AUC scores vector from cells_AUC_cal
auc_matrix <- assay(cells_AUC_cal)
gc_annot$OXPHOS_AUC <- auc_matrix["OXPHOS", colnames(gc_annot)]



#Classify cells into High and Low OXPHOS groups
gc_annot$OXPHOS_Group <- ifelse(gc_annot$OXPHOS_AUC > 0.55, "High", "Low")
gc_annot$OXPHOS_Group <- factor(gc_annot$OXPHOS_Group, levels = c("Low", "High"))



#umap plot of oxphos high and low
DimPlot(gc_annot, group.by = "OXPHOS_Group", reduction = "umap", pt.size = 0.5, cols = c("gray80", "firebrick")) +
  ggtitle("UMAP: OXPHOS High vs Low Cells") +
  theme(plot.title = element_text(hjust = 0.5))




ggplot(gc_annot@meta.data, aes(x = OXPHOS_AUC, fill = OXPHOS_Group)) +
  geom_histogram(bins = 40, alpha = 0.7, position = "identity") +
  scale_fill_manual(values = c("Low" = "blue", "High" = "firebrick")) +
  theme_minimal() +
  labs(title = "Histogram of OXPHOS AUCell Scores", x = "OXPHOS AUCell Score", y = "Cell Count")





#save the plot
ggsave('plots/Hist_oxphos_AUCell.png', width = 8, height = 6)

#plot a histogram of AUC scores
hist(gc_annot$Osphos_score, breaks = 70, col = "skyblue", main = "AUC Score Distribution", xlab = "MRG AUC Score")
abline(v = 0.5, col = "red", lty = 2)


#GO and KEGG enrichment analysis for DEGs between High vs Low OXPHOS groups.
library(clusterProfiler)
library(org.Hs.eg.db)     # For human gene annotation
library(enrichplot)       # For GO/KEGG visualization

#set cell identities and find DEGs
Idents(gc_annot) <- "OS_High"

# Find DEGs between High and Low OXPHOS score cells
degs <- FindMarkers(gc_annot, ident.1 = "High", ident.2 = "Low", assay = "SCT", slot = "data", logfc.threshold = 0.25)

# Optional: View top DEGs
head(degs)

#get significant genes
sig_genes <- degs %>%
  filter(p_val_adj < 0.05) %>%
  rownames_to_column(var = "gene")
view(sig_genes)

gene_entrez <- bitr(sig_genes$gene, fromType = "SYMBOL",
                    toType = "ENTREZID", OrgDb = org.Hs.eg.db)

#perform GO enrichment analysis
ego <- enrichGO(gene          = gene_entrez$ENTREZID,
                OrgDb         = org.Hs.eg.db,
                keyType       = "ENTREZID",
                ont           = "BP",   # Biological Process
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05)

#kegg pathway enrichment analysis
ekegg <- enrichKEGG(gene         = gene_entrez$ENTREZID,
                    organism     = 'hsa',
                    pvalueCutoff = 0.05)

# Top 10 GO terms
barplot(ego, showCategory = 10, 
        title = "GO Biological Process",
        color = c("darkred", "purple", "darkblue"))


library(ggplot2)

p <- barplot(ego, showCategory = 10, title = "GO Biological Process", color = "p.adjust")

# Customize colors with scale_fill_gradient (or scale_fill_gradientn for multiple colors)
p + scale_fill_gradient(low = "red", high = "darkblue")

#ssave the plot
ggsave('plots/GO_bp.png', width = 8, height = 6)

# KEGG
q <- barplot(ekegg, showCategory = 10, title = "KEGG Pathways")
q + scale_fill_gradient(low = "red", high = "darkblue")

















###################################################################3
#read ferroptosis suppressor genes
ferroptosis_suppressor <- read_csv("ferroptosis_suppressor.csv")
View(ferroptosis_suppressor)

#filter Human genes from ferroptosis suppressor
ferroptosis_suppressor_human <- ferroptosis_suppressor %>%
  filter(testin == "Human")
view(ferroptosis_suppressor_human)



#now intersect
matching_genes_ferroptosis <- intersect(all_markers$gene, ferroptosis_suppressor_human$symbol)
view(matching_genes_ferroptosis)

matching_genes_ferroptosis <- as.data.frame(matching_genes_ferroptosis)
colnames(matching_genes_ferroptosis) <- "Genes"

matching_markers_ferroptosis <- all_markers[all_markers$gene %in% matching_genes_ferroptosis$Genes, ]
View(matching_markers_ferroptosis)






















# If identities are annotated
Idents(gc_annot) <- "Broad_Annotation"
epi <- subset(gc_annot, idents = "Epithelial")






#recluster only the epithelial cells
epi <- NormalizeData(epi)
epi <- FindVariableFeatures(epi)
epi <- ScaleData(epi)
epi <- RunPCA(epi)
epi <- FindNeighbors(epi, dims = 1:20)
epi <- FindClusters(epi, resolution = 0.15)   # adjust resolution if needed
epi <- RunUMAP(epi, dims = 1:20)
DimPlot(epi, label = TRUE)


matching_genes[!matching_genes %in% rownames(epi)]  # shows missing genes
library(patchwork)



# Plot the matching genes on the UMAP
FeaturePlot(epi, features = matching_genes$Genes[17:22], cols = c("lightgrey", "blue"), ncol = 3)
view(matching_genes)




#run tsne plot
DimPlot(epi, reduction = "tsne", label = TRUE, pt.size = 0.5) +
  ggtitle("t-SNE Plot of Epithelial Cells")

#number of cells in eacg cluster
table(epi$seurat_clusters)


#save updated seurat
saveRDS(epi, file = "RDS/epi_reclustered.rds")

Idents(epi) <- "seurat_clusters"  # Set identity to new epi clusters


# Calculate average OXPHOS score
epi@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(mean_oxphos = mean(OXPHOS_AUC)) %>%
  arrange(desc(mean_oxphos))
VlnPlot(epi, features = "OXPHOS_AUC", group.by = "seurat_clusters", pt.size = 0.1)

view(epi@meta.data)
hist(epi$OXPHOS_AUC, breaks = 50, main = "Oxphos Score Distribution", xlab = "oxphos Score")


#use quantile cutoff as the data is left skewed
q75 <- quantile(epi$OXPHOS_AUC, 0.75)
q25 <- quantile(epi$OXPHOS_AUC, 0.25)

epi$Oxphos_Group <- ifelse(epi$OXPHOS_AUC >= q75, "High",
                           ifelse(epi$OXPHOS_AUC <= q25, "Low", "Intermediate"))


# Add the Oxphos_Group metadata to the Seurat object metadata slot
epi@meta.data$Oxphos_Group <- epi$Oxphos_Group

# Plot UMAP colored by Oxphos_Group
DimPlot(epi, group.by = "Oxphos_Group", pt.size = 0.5) +
  scale_color_manual(values = c("Low" = "grey", "Intermediate" = "grey", "High" = "red")) +
  ggtitle("UMAP Plot colored by OXPHOS Activity Groups")
#save the plot
ggsave('plots/umap_epi_oxphos_groups_dis.png', width = 8, height = 6)


#plot histogram
hist(epi$OXPHOS_AUC, breaks = 50, main = "Oxphos Score Distribution", xlab = "oxphos Score")





#subset and recluster high OXPHOS epithelial cells to see if they form functional subgroups:
epi_high <- subset(epi, subset = OXPHOS_AUC > median(epi$OXPHOS_AUC))
epi_high <- ScaleData(epi_high)
epi_high <- RunPCA(epi_high)
epi_high <- FindNeighbors(epi_high, dims = 1:20)
epi_high <- FindClusters(epi_high, resolution = 0.2)
epi_high <- RunUMAP(epi_high, dims = 1:20)

DimPlot(epi_high, label = TRUE) + ggtitle("Clusters Within High OXPHOS Epithelial Cells")
table(epi_high$seurat_clusters)


#save with ggsave
ggsave('plots/clusters_within_high_oxphos_epi.png', width = 8, height = 6)

#analyze functional state for high oxphos epithelial clusters
#find cluster specific markers
epi_high <- PrepSCTFindMarkers(epi_high)
markers_epi_high <- FindAllMarkers(
  epi_high,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)









#run pathway analysis
library(clusterProfiler)
library(org.Hs.eg.db)

# Example for cluster 0:
cluster0_genes <- markers_epi_high %>%
  dplyr::filter(cluster == 0) %>%
  dplyr::pull(gene)

ego_cluster0 <- enrichGO(
  gene = cluster0_genes,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",    # Biological Process
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)
#top 10 GO terms 
barplot(ego_cluster0, showCategory = 9 , title = "GO Biological Process for Cluster 0") +
  scale_fill_gradient(low = "red", high = "darkblue")
# View the results
dotplot(ego_cluster0, showCategory = 10)

#do it for cluster 1
cluster1_genes <- markers_epi_high %>%
  dplyr::filter(cluster == 1) %>%
  dplyr::pull(gene)
ego_cluster1 <- enrichGO(
  gene = cluster1_genes,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",    # Biological Process
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)
#barplot for cluster 1  
barplot(ego_cluster1, showCategory = 9 , title = "GO Biological Process for Cluster 1") +
  scale_fill_gradient(low = "red", high = "darkblue")






#save the plot
ggsave('plots/GO_bp_cluster1_high_oxphos.png', width = 8, height = 6)


#checking cell cycle state
# Load built-in gene lists
DefaultAssay(epi_high) <- "SCT"   # or "SCT" if you used SCTransform

# Normalize if needed (skip if already normalized)
epi_high <- SCTransform(epi_high)
epi_high <- FindVariableFeatures(epi_high)

# Run cell cycle scoring again with genes present in your object
s.genes.present <- s.genes[s.genes %in% rownames(epi_high)]
g2m.genes.present <- g2m.genes[g2m.genes %in% rownames(epi_high)]

epi_high <- CellCycleScoring(
  object = epi_high,
  s.features = s.genes.present,
  g2m.features = g2m.genes.present,
  set.ident = TRUE
)



# On UMAP
DimPlot(epi_high, group.by = "Phase", label = TRUE) + ggtitle("Cell Cycle Phase")
table(epi_high$seurat_clusters, epi_high$Phase)

#save with ggsave
ggsave('plots/umap_epi_high_cell_cycle_phase.png', width = 8, height = 6)



#Get curated gene sets
# Install if you don't have it

library(GSEABase)

gmt_file <- "h.all.v2025.1.Hs.symbols.gmt"  # replace with your file path

gene_sets <- getGmt(gmt_file)



#extract genes for each pathway of interest
stress_genes <- unique(c(
  geneIds(gene_sets[["HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY"]]),
  geneIds(gene_sets[["HALLMARK_UNFOLDED_PROTEIN_RESPONSE"]])
))

inflammation_genes <- unique(c(
  geneIds(gene_sets[["HALLMARK_TNFA_SIGNALING_VIA_NFKB"]]),
  geneIds(gene_sets[["HALLMARK_INFLAMMATORY_RESPONSE"]]),
  geneIds(gene_sets[["HALLMARK_IL6_JAK_STAT3_SIGNALING"]])
))

apoptosis_genes <- geneIds(gene_sets[["HALLMARK_APOPTOSIS"]])



#keep only genes present in my dataset
stress_genes <- intersect(stress_genes, rownames(epi_high))
inflammation_genes <- intersect(inflammation_genes, rownames(epi_high))
apoptosis_genes <- intersect(apoptosis_genes, rownames(epi_high))



#add module score to epi_high
epi_high <- AddModuleScore(epi_high, features = list(stress_genes), name = "Stress_Score")
epi_high <- AddModuleScore(epi_high, features = list(inflammation_genes), name = "Inflammation_Score")
epi_high <- AddModuleScore(epi_high, features = list(apoptosis_genes), name = "Apoptosis_Score")


#visualize score
FeaturePlot(epi_high, features = c("Stress_Score1", "Inflammation_Score1", "Apoptosis_Score1"), cols = c("lightgrey", "red"))
VlnPlot(epi_high, features = c("Stress_Score1", "Inflammation_Score1", "Apoptosis_Score1"), group.by = "seurat_clusters")


#save the scores to metadata
epi_high$Stress_Score <- epi_high$Stress_Score1
epi_high$Inflammation_Score <- epi_high$Inflammation_Score1
epi_high$Apoptosis_Score <- epi_high$Apoptosis_Score1




#check histogram for distribution for choosing threshold
hist(epi_high$Stress_Score1, breaks = 50, main = "Stress Score Distribution", xlab = "Stress Score")
hist(epi_high$Inflammation_Score1, breaks = 50, main = "Inflammation Score Distribution", xlab = "Inflammation Score")
hist(epi_high$Apoptosis_Score1, breaks = 50, main = "Apoptosis Score Distribution", xlab = "Apoptosis Score")



#define threshold
# Calculate 25th and 75th percentiles for Stress_Score1
q75 <- quantile(epi_high$Stress_Score1, 0.75)
q25 <- quantile(epi_high$Stress_Score1, 0.25)

# Assign Stress Groups based on quantiles
epi_high$Stress_Group <- ifelse(epi_high$Stress_Score1 >= q75, "High",
                                ifelse(epi_high$Stress_Score1 <= q25, "Low", "Intermediate"))


table(epi_high$Stress_Group)
hist(epi_high$Stress_Score1, breaks = 50, main = "Stress Score Distribution", xlab = "Stress Score")
abline(v = mean_score, col = "blue", lwd = 2)
abline(v = mean_score + sd_score, col = "red", lwd = 2, lty = 2)
abline(v = mean_score - sd_score, col = "red", lwd = 2, lty = 2)



# Stress example
DimPlot(epi_high, group.by = "Stress_Group") +
  scale_color_manual(values = c("Low" = "grey", "Intermediate" = "grey", "High" = "firebrick")) +
  ggtitle("Stress Score: High vs Low")



# Calculate 25th and 75th percentiles for Stress_Score1
q75 <- quantile(epi_high$Inflammation_Score1, 0.75)
q25 <- quantile(epi_high$Inflammation_Score1, 0.25)

# Assign Stress Groups based on quantiles
epi_high$Inflammation_Group <- ifelse(epi_high$Inflammation_Score1 >= q75, "High",
                                ifelse(epi_high$Inflammation_Score1 <= q25, "Low", "Intermediate"))


DimPlot(epi_high, group.by = "Inflammation_Group") +
  scale_color_manual(values = c("Low" = "grey", "Intermediate" = "grey", "High" = "firebrick")) +
  ggtitle("Inflammation Score: High vs Low")


# Calculate 25th and 75th percentiles for Stress_Score1
q75 <- quantile(epi_high$Apoptosis_Score1, 0.75)
q25 <- quantile(epi_high$Apoptosis_Score1, 0.25)

# Assign Stress Groups based on quantiles
epi_high$Apoptosis_Group <- ifelse(epi_high$Apoptosis_Score1 >= q75, "High",
                                      ifelse(epi_high$Apoptosis_Score1 <= q25, "Low", "Intermediate"))


DimPlot(epi_high, group.by = "Apoptosis_Group") +
  scale_color_manual(values = c("Low" = "grey", "Intermediate" = "grey", "High" = "firebrick")) +
  ggtitle("Apoptosis Score: High vs Low")

#save with ggsave
ggsave('plots/Inflammation_score_high.png', width = 8, height = 6)
ggsave('plots/Apoptosis_score.png', width = 8, height = 6)




