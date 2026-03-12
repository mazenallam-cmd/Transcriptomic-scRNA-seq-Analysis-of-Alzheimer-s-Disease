library(Seurat)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggrepel)


#read_data

counts <- readRDS(file.choose())
metadata <- readRDS(file.choose())
sample_metadata <- readRDS(file.choose())
gene_metadata <- readRDS(file.choose())





#create_seurat_object
cell_full_info <- merge(metadata, sample_metadata, by = "sample_id")
rownames(cell_full_info) <- cell_full_info$cell_id
cell_full_info <- cell_full_info[colnames(counts), ]

sn_obj <- CreateSeuratObject(counts = counts, meta.data = cell_full_info, project = "SingleCellAnalysis")


sn_obj

# QC
range(sn_obj$percent.mt)
range(sn_obj$nCount_RNA)
range(sn_obj$nFeature_RNA)

VlnPlot(sn_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


Fsn_obj <- subset(sn_obj, subset = nCount_RNA > 500 & 
                    nCount_RNA < 6000 & 
                    percent.mt < 5)
range(Fsn_obj$nFeature_RNA)



#normalization
Fsn_obj <- NormalizeData(Fsn_obj, normalization.method = "LogNormalize", scale.factor = 10000)
Fsn_obj
#variable features
Fsn_obj <- FindVariableFeatures(Fsn_obj, selection.method = "vst", nfeatures = 2000)
top_10 <- head(VariableFeatures(Fsn_obj), 10)

#variable_plot
plot_1<-VariableFeaturePlot(Fsn_obj)
plot_2<-LabelPoints(plot = plot_1, points = top_10, repel = TRUE)
plot_2

Fsn_obj <- ScaleData(Fsn_obj, features = VariableFeatures(Fsn_obj), 
                     vars.to.regress = c( "percent.mt"))
Fsn_obj <- RunPCA(Fsn_obj, features = VariableFeatures(Fsn_obj))
ElbowPlot(Fsn_obj)



#neighbor graph and cluster
Fsn_obj <- FindNeighbors(Fsn_obj, dims = 1:20)
Fsn_obj <- FindClusters(Fsn_obj, resolution = 0.5)

Fsn_obj <- RunUMAP(Fsn_obj, dims = 1:20)

FeaturePlot(Fsn_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.5)
view(Fsn_obj)



# Data frame: cell barcode, cluster, cell_type
df <- data.frame(
  cell_id = Fsn_obj$cell_id,
  cluster = Idents(Fsn_obj),
  cell_type = Fsn_obj$cell_type
)

# Find mode per cluster
mode_cell_type <- df %>% group_by(cluster) %>% 
  count(cell_type) %>%
  top_n(1, n) %>%           
  slice(1) %>%                
  ungroup()

# Create a named vector for renaming
cluster_names <- setNames(mode_cell_type$cell_type, mode_cell_type$cluster)

Fsn_obj <- RenameIdents(Fsn_obj, cluster_names)
DimPlot(Fsn_obj, reduction = "umap", label = TRUE, pt.size = 0.5)
DimPlot(Fsn_obj, reduction = "pca", label = TRUE, pt.size = 0.5)

# find_top10_markers_per_cluster

all_top10_markers <- list()

for (cluster in cluster_names) {
  markers <- FindMarkers(Fsn_obj, ident.1 = cluster, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  
  top10 <- markers %>% 
    rownames_to_column("gene") %>%
    top_n(n = 10, wt = avg_log2FC)
  
  all_top10_markers[[cluster]] <- top10
  
  # Filter genes that exist in the object
  valid_genes <- top10$gene[top10$gene %in% rownames(Fsn_obj)]
  
  if (length(valid_genes) == 0) {
    cat("No valid genes found for cluster:", cluster, "\n")
    next
  }
  
  # Heatmap
  p1 <- DoHeatmap(Fsn_obj, features = valid_genes) +
    ggtitle(paste("Top 10 Markers -", cluster))
  p1
  
  # Save Heatmap as PDF
  pdf(paste0("Heatmap_", cluster, ".pdf"), width = 10, height = 8)
  p1
  dev.off()
  
  # DotPlot
  p2 <- DotPlot(Fsn_obj, features = valid_genes) +
    RotatedAxis() +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank()
    ) +
    ggtitle(paste("Top 10 Markers DotPlot -", cluster))
  print(p2)
  
  # Save DotPlot as PDF
  pdf(paste0("DotPlot_", cluster, ".pdf"), width = 10, height = 8)
  print(p2)
  dev.off()
  
  # Volcano Plot
  volcano_data <- markers %>% 
    rownames_to_column("gene") %>%
    mutate(is_sig = p_val < 0.05 & avg_log2FC > 0.25)
  
  p3 <- ggplot(volcano_data, aes(x = avg_log2FC, y = -log10(p_val))) +
    geom_point(aes(color = is_sig), size = 2, alpha = 0.6) +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "gray"),
                       labels = c("TRUE" = "Significant", "FALSE" = "Not Sig")) +
    geom_vline(xintercept = 0.25, linetype = "dashed", color = "blue", alpha = 0.5) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue", alpha = 0.5) +
    geom_text_repel(data = subset(volcano_data, is_sig),
                    aes(label = gene), size = 3, max.overlaps = 20,
                    box.padding = 0.5, point.padding = 0.5) +
    labs(title = paste("Volcano Plot -", cluster),
         x = "Log2 Fold Change", y = "-Log10(p-value)", color = "Significance") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))
  print(p3)
  
  # Save Volcano Plot as PDF
  pdf(paste0("VolcanoPlot_", cluster, ".pdf"), width = 10, height = 8)
  print(p3)
  dev.off()
}

# Combine all top 10 markers into one data frame
top10_all_clusters <- do.call(rbind, all_top10_markers)
rownames(top10_all_clusters) <- NULL

# Save to CSV
write.csv(top10_all_clusters, "top10_markers_all_clusters.csv", row.names = FALSE)




table(Fsn_obj$cell_type, Fsn_obj$group)

table(Fsn_obj$cell_type, Fsn_obj$sex.x)

table(Fsn_obj$group, Fsn_obj$sex.x)


DimPlot(Fsn_obj, reduction = "umap",label = TRUE, split.by = "sex.x", pt.size = 0.5)

DimPlot(Fsn_obj, reduction = "umap",label = TRUE, split.by = "group", pt.size = 0.5)



#DEGs for each cell type between AD and CTRL

Idents(Fsn_obj) <- Fsn_obj$cell_type
cell_types <- unique(Idents(Fsn_obj))

de_results <- list()

for (cell_type in cell_types) {
  cat("\nAnalyzing cell type:", cell_type)
  
  cells_subset <- subset(Fsn_obj, idents = cell_type)
  groups_present <- table(cells_subset$group)
  
  tryCatch({
    de_genes <- FindMarkers(
      cells_subset,
      ident.1 = "AD",
      ident.2 = "CTRL",
      group.by = "group",
      test.use = "wilcox",
      min.pct = 0.05
    )
    de_genes$gene <- rownames(de_genes)
    de_genes$cell_type <- cell_type
    de_results[[cell_type]] <- de_genes
    
  }, error = function(e) cat("  Skipping", cell_type, ":", e$message, "\n"))
}

de_all <- bind_rows(de_results)
table(de_all$cell_type)

de_all <- de_all %>%
  mutate(sig_threshold = p_val < 0.05)

top_de_genes <- de_all %>%
  group_by(cell_type) %>%
  arrange(p_val) %>%
  slice_head(n = 20)
print(top_de_genes)

write.csv(de_all, "DEG_per_all_celltypes_AD_vs_CTRL.csv")


# VOLCANO PLOTS for each cell type

pdf("Volcano_Plots_AD_vs_CTRL.pdf", width = 12, height = 8)

for (cell_type in names(de_results)) {
  de_data <- de_results[[cell_type]]
  de_data$is_sig <- de_data$p_val < 0.05 & abs(de_data$avg_log2FC) > 0.25
  
  p <- ggplot(de_data, aes(x = avg_log2FC, y = -log10(p_val))) +
    geom_point(aes(color = is_sig), size = 2, alpha = 0.6) +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "gray"),
                       labels = c("TRUE" = "Significant", "FALSE" = "Not Sig")) +
    geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "blue", alpha = 0.5) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue", alpha = 0.5) +
    geom_text_repel(data = subset(de_data, is_sig),
                    aes(label = gene), size = 3, max.overlaps = 20,
                    box.padding = 0.5, point.padding = 0.5) +
    labs(title = paste("Volcano Plot -", cell_type, "(AD vs CTRL)"),
         x = "Log2 Fold Change", y = "-Log10(p-value)", color = "Significance") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))
  print(p)
}

dev.off() 

# HEATMAP of top DE genes for each cell type
pdf("Heatmaps_DE_genes_AD_vs_CTRL.pdf", width = 12, height = 10)

for (cell_type in names(de_results)) {
  cat("  Processing heatmap for:", cell_type, "\n")
  
  top_genes <- de_results[[cell_type]] %>%
    arrange(p_val) %>%
    slice_head(n = 30) %>%
    pull(gene)
  cells_subset <- subset(Fsn_obj, idents = cell_type)
  
  p <- DoHeatmap(cells_subset, features = top_genes, group.by = "group",
                 label = TRUE, size = 3) +
    ggtitle(paste("Top 30 DE Genes -", cell_type)) +
    theme(title = element_text(size = 12, face = "bold"))
  print(p)
}

dev.off()  


# AGE_ANALYSIS
Fsn_obj$age_interval <- cut(Fsn_obj$age,
                            breaks = c(59, 70, 80, 90),
                            labels = c("60-70", "71-80", "81-90"),
                            include.lowest = TRUE)

table(Fsn_obj$age_interval, Fsn_obj$group)

Idents(Fsn_obj) <- Fsn_obj$cell_type

age_intervals <- c("60-70", "71-80", "81-90")
de_age_results <- list()

for (interval in age_intervals) {
  
  obj_age <- subset(Fsn_obj, subset = age_interval == interval)
  print(table(obj_age$age_interval, obj_age$group))  
  
  de_interval_results <- list()
  
  for (cell_type in cell_types) {  
    
    cells_subset <- subset(obj_age, idents = cell_type)
    groups_present <- table(cells_subset$group)
    
    tryCatch({  
      de_genes <- FindMarkers(
        cells_subset,
        ident.1 = "AD",
        ident.2 = "CTRL",
        group.by = "group",
        test.use = "wilcox",
        min.pct = 0.05
      )
      de_genes$gene <- rownames(de_genes)
      de_genes$cell_type <- cell_type
      de_genes$age_interval <- interval
      de_interval_results[[cell_type]] <- de_genes  
      
    }, error = function(e) cat("  Skipping", cell_type, "in", interval, ":", e$message, "\n")) 
  }  
  
  de_age_results[[interval]] <- bind_rows(de_interval_results)  
}  

de_age_all <- bind_rows(de_age_results)
table(de_age_all$age_interval, de_age_all$cell_type)

write.csv(de_age_all, "DEGs_Age_AD_vs_CTRL.csv" )


top_n_genes_age <- 30

for (interval in names(de_age_results)) {
  de_interval <- de_age_results[[interval]]
  
  pdf(paste0("Heatmaps_DE_genes_Age_", gsub("[^A-Za-z0-9]", "_", interval), ".pdf"), width = 12, height = 10)
  
  celltypes_in_interval <- unique(de_interval$cell_type)
  for (cell_type in celltypes_in_interval) {
    de_data_ct <- de_interval %>% filter(cell_type == !!cell_type)
    
    if (nrow(de_data_ct) == 0) next
    
    # Select top genes by p-value
    top_genes <- de_data_ct %>%
      arrange(p_val) %>%
      slice_head(n = top_n_genes_age) %>%
      pull(gene)
    
   
    valid_genes <- top_genes[top_genes %in% rownames(Fsn_obj)]
    cells_subset <- subset(Fsn_obj, subset = age_interval == interval)
    cells_subset <- subset(cells_subset, idents = cell_type)
   
    
    p <- DoHeatmap(cells_subset, features = valid_genes, group.by = "group",
                   label = TRUE, size = 3) +
      ggtitle(paste("Top", top_n_genes_age, "DE Genes -", cell_type, "-", interval)) +
      theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))
    
    print(p)
  }
  
  dev.off()
}




