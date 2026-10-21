# 1. Introductory Comments:
# Shamim_Roshni - Full Pipeline for 12 rds files
# Clearly describes the purpose of the code: processing and quality control (QC) for 12 single-cell RNA-seq .rds files.
# This code worked perfectly and saves RDS data files after QC
# This code is good to save RDS files after QC; plots generated here are not our desired plots
# Forget about these plots—they are not our desired final plots, just for my test

# 2. Memory Limit Setting:
# Increases the global memory allocation limit to handle large datasets (400 GB specified here).
options(future.globals.maxSize = 10 * 1024^3)  # Set to 400GB (increase if needed)

# 3. Package Loading and Installation:
# Load required libraries; install them if not present
if (!requireNamespace("Seurat", quietly = TRUE)) install.packages("Seurat")
if (!requireNamespace("SeuratObject", quietly = TRUE)) install.packages("SeuratObject")
if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) BiocManager::install("SingleCellExperiment")
if (!requireNamespace("scDblFinder", quietly = TRUE)) BiocManager::install("scDblFinder")

# Checks for the required libraries (Seurat, ggplot2, dplyr, etc.) and installs them automatically if missing.
# Loads required libraries, suppressing unnecessary startup messages.
library(Seurat)
library(SeuratObject)
library(patchwork)
library(ggplot2)
library(SingleCellExperiment)
library(scDblFinder)

# 4. Directory Setup:
# Defines data and QC results folders; creates QC_Results4 if it doesn't exist.

# For mac path of all data
#data_folder <- "/Users/abmshamimulhasan/"

# For RCC/midway3 path
data_folder <- "/home/abmshamimuh/"
qc_folder <- file.path(data_folder, "QC_Results_RDS_DATA")
if (!dir.exists(qc_folder)) dir.create(qc_folder)

# 5. Defining Sample Names:
# Lists the 12 samples to be processed.
samples <- c(
  "GSM4339769_C141", "GSM4339770_C142", "GSM4339771_C143",
  "GSM4339772_C144", "GSM4339773_C145", "GSM4339774_C146",
  "GSM4475048_C51", "GSM4475049_C52", "GSM4475050_C100",
  "GSM4475051_C148", "GSM4475052_C149", "GSM4475053_C152"
)


# 6. Function to Process Each Sample:
# Process and QC each sample - existing function
# Defined as a reusable function for QC and preprocessing for each sample individually:

# Check for Existing Processed Data:
# If a sample is already processed, skip reprocessing.
process_sample <- function(sample_name) {
  final_rds <- file.path(qc_folder, paste0(sample_name, "_seurat_final.rds"))
  
  if (file.exists(final_rds)) {
    cat("\nSkipping already processed sample:", sample_name, "\n")
    return(readRDS(final_rds))
  }
  
  cat("\nProcessing sample:", sample_name, "\n")
  
  # Reading and Creating Seurat Object:
  # Loads data and initializes the Seurat object.
  filepath <- file.path(data_folder, paste0(sample_name, "_filtered_feature_bc_matrix.h5"))
  seurat_obj <- Read10X_h5(filepath)
  seurat_obj <- CreateSeuratObject(counts = seurat_obj, project = sample_name)
  
  # Mitochondrial QC metric:Calculate Mitochondrial Gene Content:
  # Computes mitochondrial gene percentage for QC purposes.
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  
  # QC Violin plot:
  # Creates and saves QC plots for each sample.
  vln <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  ggsave(file.path(qc_folder, paste0(sample_name, "_QC_ViolinPlot.png")), plot = vln, width = 10, height = 5)
  
  # Apply QC Filters: Filter Low-Quality Cells:
  # Keeps only high-quality cells based on gene count and mitochondrial percentage.
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
  
  # Normalization & dimensional reduction:Normalizes the data using SCTransform, performs PCA and UMAP.
  seurat_obj <- SCTransform(seurat_obj, verbose = FALSE)
  seurat_obj <- RunPCA(seurat_obj)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)
  
  # Doublet Detection (using scDblFinder)
  # Identifies doublets (cells likely representing more than one cell).
  sce <- as.SingleCellExperiment(seurat_obj)
  sce <- scDblFinder(sce)
  seurat_obj$scDblFinder_class <- sce$scDblFinder.class
  
  # Doublet Visualization & Removal: Removes doublets from the dataset.
  umap_doublet <- DimPlot(seurat_obj, reduction = "umap", group.by = "scDblFinder_class") +
    ggtitle(paste0(sample_name, " - Doublet Classification (scDblFinder)"))
  ggsave(file.path(qc_folder, paste0(sample_name, "_scDblFinder_UMAP.png")), plot = umap_doublet, width = 8, height = 6)
  
  seurat_obj <- subset(seurat_obj, subset = scDblFinder_class == "singlet")
  
  # Clustering & Saving the Final Object: Clusters the cells and saves the resulting Seurat object.
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30)
  seurat_obj <- FindClusters(seurat_obj, resolution = 1.2)
  
  umap_final <- DimPlot(seurat_obj, reduction = "umap", label = TRUE) +
    ggtitle(paste0(sample_name, " - Final Clustering (After Doublet Removal)"))
  ggsave(file.path(qc_folder, paste0(sample_name, "_Final_UMAP.png")), plot = umap_final, width = 8, height = 6)
  
  # Saves the resulting Seurat object. This is our desired rds files after QC
  saveRDS(seurat_obj, final_rds)
  return(seurat_obj)
}


# Additional Part, not our desired plots

# 7. Process All Samples: Processes each sample using the above-defined function.
# Process all samples (load existing if already processed)
seurat_list <- lapply(samples, process_sample)
cat("\nAll samples processed or loaded from cache.\n")

# Add sample metadata
for (i in seq_along(seurat_list)) {
  seurat_list[[i]]$sample <- samples[i]
}

# 8. Data Integration (Across all Samples):
# Integrates individual Seurat objects into a single, combined dataset.
features <- SelectIntegrationFeatures(object.list = seurat_list, nfeatures = 3000)
seurat_list <- PrepSCTIntegration(object.list = seurat_list, anchor.features = features)

anchors <- FindIntegrationAnchors(object.list = seurat_list, normalization.method = "SCT", anchor.features = features)
combined_seurat <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

# Standard processing on combined object
combined_seurat <- RunPCA(combined_seurat)
combined_seurat <- RunUMAP(combined_seurat, dims = 1:30)

# 9. Clustering on Combined Data (target 14 clusters if possible):
target_clusters <- 14
found <- FALSE
for (res in seq(0.1, 2.0, by = 0.1)) {
  combined_seurat <- FindNeighbors(combined_seurat, dims = 1:30)
  combined_seurat <- FindClusters(combined_seurat, resolution = res)
  if (length(unique(Idents(combined_seurat))) == target_clusters) {
    cat(paste0("Found exactly ", target_clusters, " clusters at resolution = ", res, "\n"))
    found <- TRUE
    break
  }
}
if (!found) {
  cat("Could not find exactly 14 clusters. Using resolution = 1.0.\n")
  combined_seurat <- FindClusters(combined_seurat, resolution = 1.0)
}

# 10. Saving Final Integrated Object:Saves the combined, clustered dataset.
saveRDS(combined_seurat, file.path(qc_folder, "Combined_Seurat_Object.rds"))

# 11.Plots:
# Plot 1 - Combined UMAP Plot (14 clusters)
# 1. Combined UMAP
p1 <- DimPlot(combined_seurat, reduction = "umap", label = TRUE, pt.size = 0.5) +
  ggtitle("Combined UMAP - All 12 Samples (14 Clusters)")
ggsave(file.path(qc_folder, "Combined_12_Samples_UMAP.png"), plot = p1, width = 10, height = 6)

# Plot 2 - Sample Representation Per Cluster
# 2. Sample Proportion Per Cluster
combined_seurat$sample <- sapply(strsplit(rownames(combined_seurat@meta.data), "_"), `[`, 1)

# Count sample contributions per cluster
cluster_sample_table <- table(Idents(combined_seurat), combined_seurat$sample)
cluster_sample_df <- as.data.frame(cluster_sample_table)
colnames(cluster_sample_df) <- c("Cluster", "Sample", "Count")

p2 <- ggplot(cluster_sample_df, aes(x = Cluster, y = Count, fill = Sample)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(y = "Proportion", title = "Proportion of Each Sample in Each Cluster") +
  theme_minimal() +
  coord_flip()

ggsave(file.path(qc_folder, "Sample_Proportion_Per_Cluster.png"), plot = p2, width = 10, height = 6)

# Plot 3 - GZMK Expression Per Cluster (Highest to Lowest)
# Calculate average expression of GZMK in each cluster
# 3. GZMK Expression Per Cluster
avg_exp <- AverageExpression(combined_seurat, features = "GZMK", group.by = "seurat_clusters")
avg_exp_df <- as.data.frame(avg_exp$SCT)
avg_exp_df$Cluster <- rownames(avg_exp_df)

p3 <- ggplot(avg_exp_df, aes(x = reorder(Cluster, -GZMK), y = GZMK)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title = "GZMK Expression Across Clusters", x = "Cluster", y = "Average GZMK Expression") +
  theme_minimal()

ggsave(file.path(qc_folder, "GZMK_Expression_Per_Cluster.png"), plot = p3, width = 10, height = 6)
