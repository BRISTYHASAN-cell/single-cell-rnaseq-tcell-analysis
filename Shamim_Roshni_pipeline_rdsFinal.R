# 1. Introductory Comments:
# Shamim_Roshni - Full Pipeline for 12 rds files
# Clearly describes the purpose of the code: processing and quality control (QC) for 12 single-cell RNA-seq .rds files.
# This code worked perfectly and saves RDS data files after QC
# This code is good to save RDS files after QC; plots generated here are not our desired plots
# Forget about these plots—they are not our desired final plots, just for my test

# 2. Memory Limit Setting:
# Increases the global memory allocation limit to handle large datasets (400 GB specified here).
options(future.globals.maxSize = 8 * 1024^3)  # Set to 300 GB to handle larger samples

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
# Create QC output directory
qc_folder <- file.path(data_folder, "QC_Results_12RDS")
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
  
  # Skip already processed samples
  if (file.exists(final_rds)) {
    cat("\n kipping already processed sample:", sample_name, "\n")
    return(readRDS(final_rds))
  }
  
  cat("\n Processing sample:", sample_name, "\n")
  
  # Loading 10X Single-cell Data
  # Reading and Creating Seurat Object:
  # Loads data and initializes the Seurat object.
  filepath <- file.path(data_folder, paste0(sample_name, "_filtered_feature_bc_matrix.h5"))
  seurat_obj <- Read10X_h5(filepath)
  seurat_obj <- CreateSeuratObject(counts = seurat_obj, project = sample_name)
  
  # Mitochondrial QC metric:Calculate Mitochondrial Gene Content:
  # Calculating Percentage of Mitochondrial Genes
  # Computes the mitochondrial gene expression percentage for quality control purposes.
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  
  # QC Violin plot:Generates and saves violin plots for key QC metrics (nFeature_RNA, nCount_RNA, and mitochondrial percentage).
  # Creates and saves QC plots for each sample.
  vln <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  ggsave(file.path(qc_folder, paste0(sample_name, "_QC_ViolinPlot.png")), plot = vln, width = 10, height = 5)
  
  # Apply QC filters (we can modify thresholds if needed)
  # Apply QC Filters: Filter Low-Quality Cells:Filters out cells with too few/many genes or high mitochondrial content.
  # Keeps only high-quality cells based on gene count and mitochondrial percentage.
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
  
  # Preprocessing - Normalization and Dimensionality Reduction
  # Performs normalization (SCTransform), PCA for dimensionality reduction, and UMAP visualization.
  seurat_obj <- SCTransform(seurat_obj, verbose = FALSE)
  seurat_obj <- RunPCA(seurat_obj)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)
  
  # Identifies doublets (cells likely representing more than one cell).
  # Converts the data into a SingleCellExperiment object for doublet detection using scDblFinder.
  sce <- as.SingleCellExperiment(seurat_obj)
  sce <- scDblFinder(sce)
  
  # Stores doublet classification results back in the Seurat object.
  seurat_obj$scDblFinder_class <- sce$scDblFinder.class
  
  # Save Doublet Classification UMAP Plot
  # Creates and saves UMAP plots visualizing detected doublets vs. singlets.
  umap_doublet <- DimPlot(seurat_obj, reduction = "umap", group.by = "scDblFinder_class", pt.size = 0.5) +
    ggtitle(paste0(sample_name, " - Doublet Classification (scDblFinder)"))
  ggsave(file.path(qc_folder, paste0(sample_name, "_scDblFinder_UMAP.png")), plot = umap_doublet, width = 8, height = 6)
  
  # Filtering Out Detected Doublets
  # Keeps only single cells identified by the doublet-detection algorithm.
  seurat_obj <- subset(seurat_obj, subset = scDblFinder_class == "singlet")
  
  # Clustering - Try to force exactly 14 clusters if possible
  # Iteratively attempts to cluster cells into exactly 14 clusters by adjusting resolution.
  target_clusters <- 14
  found <- FALSE
  for (res in seq(0.1, 2.0, by = 0.1)) {
    seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30)
    seurat_obj <- FindClusters(seurat_obj, resolution = res)
    num_clusters <- length(unique(Idents(seurat_obj)))
    if (num_clusters == target_clusters) {
      cat(paste0("ound exactly ", target_clusters, " clusters at resolution = ", res, "\n"))
      found <- TRUE
      break
    }
  }
  
  if (!found) {
    cat("️ Could not find exactly 14 clusters. Proceeding with resolution = 1.2 (default).\n")
    seurat_obj <- FindClusters(seurat_obj, resolution = 1.2)
  }
  
  # Final UMAP plot after clustering
  # Saves the final clustering visualization.
  umap_final <- DimPlot(seurat_obj, reduction = "umap", label = TRUE, pt.size = 0.5) +
    ggtitle(paste0(sample_name, " - Final Clustering (After Doublet Removal)"))
  ggsave(file.path(qc_folder, paste0(sample_name, "_Final_UMAP.png")), plot = umap_final, width = 8, height = 6)
  
  # Identifying Marker Genes: Finds marker genes for clusters and saves them as CSV files.
  # Differential Expression (Marker genes for all clusters)
  cluster_markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  write.csv(cluster_markers, file = file.path(qc_folder, paste0(sample_name, "_ClusterMarkers.csv")), row.names = FALSE)
  
  # Save final Seurat object: Stores the processed data as an .rds file for future loading.
  # This is our desired rds files after QC
  saveRDS(seurat_obj, final_rds)
  
  return(seurat_obj)
}

# Process all samples
# Applies the entire QC and processing workflow across all defined samples.
seurat_list <- lapply(samples, process_sample)


# Final confirmation message indicating workflow completion.
cat("\n All samples processed or skipped if already completed.\n")
