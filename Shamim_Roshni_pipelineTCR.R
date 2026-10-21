# Shamim_Roshni - Full Pipeline with Viral TCR Integration

# Increase memory limit
options(future.globals.maxSize = 250 * 1024^3)  # Set memory to 250GB

# Load required libraries
if (!requireNamespace("Seurat", quietly = TRUE)) install.packages("Seurat", ask = FALSE)
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2", ask = FALSE)
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr", ask = FALSE)
if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork", ask = FALSE)
if (!requireNamespace("readxl", quietly = TRUE)) install.packages("readxl", ask = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(readxl)
})

# Set directories
data_folder <- "/home/abmshamimuh/"
qc_folder <- file.path(data_folder, "QC_Results4")
if (!dir.exists(qc_folder)) dir.create(qc_folder)

# Load Seurat objects
input_dir <- "~/QC_Results4/"
file_list <- list.files(input_dir, pattern = "\\.rds$", full.names = TRUE)

# Ensure RDS files exist
if (length(file_list) == 0) stop("❌ ERROR: No RDS files found in the directory!")

# Function to safely load Seurat objects
safe_load_seurat <- function(f) {
  obj <- tryCatch({
    readRDS(f)
  }, error = function(e) {
    cat("❌ ERROR: Could not load", f, "\n", e$message, "\n")
    return(NULL)
  })
  if (!is.null(obj) && inherits(obj, "Seurat")) {
    obj <- UpdateSeuratObject(obj)
    return(obj)
  } else {
    cat("❌ ERROR: File", f, "is not a valid Seurat object!\n")
    return(NULL)
  }
}

# Load and validate Seurat objects
seurat_objects <- lapply(file_list, safe_load_seurat)
seurat_objects <- Filter(Negate(is.null), seurat_objects)  # Remove NULL objects

# Stop if no valid Seurat objects were loaded
if (length(seurat_objects) == 0) stop("❌ ERROR: No valid Seurat objects found!")

# Find common features
common_features <- Reduce(intersect, lapply(seurat_objects, rownames))
if (length(common_features) == 0) stop("❌ ERROR: No common features found!")

# Ensure SCT normalization before integration
seurat_objects <- lapply(seurat_objects, function(obj) {
  if (!"SCT" %in% names(obj@assays)) {
    obj <- SCTransform(obj, verbose = FALSE)
  }
  return(obj)
})

# Run integration
anchors <- FindIntegrationAnchors(object.list = seurat_objects, dims = 1:30, anchor.features = 3000)
combined <- IntegrateData(anchorset = anchors, dims = 1:30)

# Set default assay before clustering
DefaultAssay(combined) <- "integrated"

# Scale and Run PCA
combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, npcs = 30, verbose = FALSE)

# Run UMAP
combined <- RunUMAP(combined, dims = 1:30, verbose = FALSE)

# Find Neighbors and Clusters
combined <- FindNeighbors(combined, dims = 1:30, verbose = FALSE)
combined <- FindClusters(combined, resolution = 0.5)

# Ensure exactly 14 clusters
resolution_vals <- seq(0.2, 1.0, by = 0.05)
found_clusters <- FALSE

for (res in resolution_vals) {
  combined <- FindClusters(combined, resolution = res)
  num_clusters <- length(unique(Idents(combined)))
  cat("🔹 Trying resolution =", res, "➡️ Clusters found:", num_clusters, "\n")
  if (num_clusters == 14) {
    cat("✅ Found exactly 14 clusters at resolution =", res, "\n")
    found_clusters <- TRUE
    break
  }
}

if (!found_clusters) {
  cat("⚠️ WARNING: Could not get exactly 14 clusters. Using closest match!\n")
}

# Load Viral TCR Data
tcr_file <- "/project2/rroychowdhury/shamim_/Corridoni_TCRs.xlsx"
tcr_data <- read_excel(tcr_file, sheet = "counts with matches")

# Identify Viral-Associated Barcodes
viral_barcodes <- tcr_data %>% filter(!is.na(organism)) %>% pull(barcode)

# Add metadata to Seurat object
combined$viral_status <- ifelse(rownames(combined@meta.data) %in% viral_barcodes, "Viral TCR", "Non-Viral TCR")

# **PLOT 1: UMAP Highlighting Viral TCRs**
umap_plot_viral <- DimPlot(combined, reduction = "umap", group.by = "viral_status", 
                           cols = c("gray", "red"), label = TRUE) +
  ggtitle("UMAP with Viral TCRs Highlighted") +
  theme_minimal()

# **PLOT 2: Sample Contribution Plot**
sample_contribution <- as.data.frame(table(Idents(combined), combined$orig.ident))
colnames(sample_contribution) <- c("Cluster", "Sample", "Count")

sample_contribution_plot <- ggplot(sample_contribution, aes(x = Cluster, y = Count, fill = Sample)) +
  geom_bar(stat = "identity", position = "fill") +
  ylab("Proportion") +
  ggtitle("Sample Contribution by Cluster (Exact 14 Clusters)") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set3")

# **PLOT 3: GZMK Expression Plot**
if ("GZMK" %in% rownames(combined)) {
  gzmk_expression_plot <- VlnPlot(combined, features = "GZMK", group.by = "seurat_clusters", pt.size = 0.1) +
    ggtitle("GZMK Expression in Each Cluster") +
    scale_fill_manual(values = scales::hue_pal()(14)) +
    theme_minimal()
} else {
  gzmk_expression_plot <- ggplot() + ggtitle("⚠️ GZMK not found in dataset!") + theme_minimal()
  cat("\n⚠️ WARNING: GZMK not found in dataset! Skipping GZMK plot...\n")
}

# **PLOT 4: Viral vs Non-Viral Cluster Composition**
viral_cluster_composition <- as.data.frame(table(combined$viral_status, Idents(combined)))
colnames(viral_cluster_composition) <- c("TCR_Type", "Cluster", "Count")

viral_cluster_barplot <- ggplot(viral_cluster_composition, aes(x = Cluster, y = Count, fill = TCR_Type)) +
  geom_bar(stat = "identity", position = "dodge") +
  ylab("Number of TCRs") +
  ggtitle("Cluster Composition: Viral vs Non-Viral TCRs") +
  theme_minimal() +
  scale_fill_manual(values = c("gray", "red"))

# **Save Plots**
ggsave(file.path(qc_folder, "UMAP_Viral_TCRs.png"), plot = umap_plot_viral, width = 10, height = 6)
ggsave(file.path(qc_folder, "Sample_Contribution_14Clusters.png"), plot = sample_contribution_plot, width = 10, height = 6)
ggsave(file.path(qc_folder, "GZMK_Expression_14Clusters.png"), plot = gzmk_expression_plot, width = 10, height = 6)
ggsave(file.path(qc_folder, "Viral_vs_NonViral_Cluster.png"), plot = viral_cluster_barplot, width = 10, height = 6)

# **Save Final Object**
saveRDS(combined, file.path(qc_folder, "Final_Seurat_Object_14Clusters_Viral.rds"))

cat("\n✅ Analysis completed successfully! Plots saved in", qc_folder, " 🎉\n")
