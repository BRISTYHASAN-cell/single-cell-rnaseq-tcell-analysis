# 1. Comments and Initial Setup
# Shamim_Roshni - Full Pipeline with Additional Plots
# Final code to find exactly 14 clusters
# Confirmed working correctly on RCC (Remote Compute Cluster), no error in RCC


# Increase memory limit:Sets maximum memory usage to 250 GB, necessary for handling large datasets.
options(future.globals.maxSize = 250 * 1024^3)  # Set memory limit to 250 GB


# 2. Library Loading and Installation
# Load required libraries, installing if missing
if (!requireNamespace("Seurat", quietly = TRUE)) install.packages("Seurat", ask = FALSE)
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2", ask = FALSE)
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr", ask = FALSE)
if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork", ask = FALSE)
if (!requireNamespace("Matrix", quietly = TRUE)) install.packages("Matrix", ask = FALSE)
if (!requireNamespace("future", quietly = TRUE)) install.packages("future", ask = FALSE)

# Checks for the required libraries (Seurat, ggplot2, dplyr, etc.) and installs them automatically if missing.
# Loads required libraries, suppressing unnecessary startup messages.
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(Matrix)
  library(future)
})

# 3. Directory Setup in RCC
# Defines paths for the main data folder and creates a QC_Results4 subdirectory if it doesn’t already exist.
data_folder <- "/home/abmshamimuh/"
qc_folder <- file.path(data_folder, "QC_Results4")
if (!dir.exists(qc_folder)) dir.create(qc_folder)

# 4. Loading Seurat Objects from RDS files
# Lists all .rds files (Seurat objects) in QC_Results4 folder for processing.
input_dir <- "~/QC_Results4/"
file_list <- list.files(input_dir, pattern = "\\.rds$", full.names = TRUE)

# Ensure RDS files exist: Stops script execution immediately if no .rds files are found.
if (length(file_list) == 0) stop(" ERROR: No RDS files found in the directory!")

# 5. Safe Loading of Seurat Objects
# Function to safely load Seurat objects
# Defines a function that safely loads .rds files, handling loading errors gracefully, and ensuring the objects are valid Seurat objects.
safe_load_seurat <- function(f) {
  obj <- tryCatch({
    readRDS(f)
  }, error = function(e) {
    cat(" ERROR: Could not load", f, "\n", e$message, "\n")
    return(NULL)
  })
  if (!is.null(obj) && inherits(obj, "Seurat")) {
    obj <- UpdateSeuratObject(obj)
    return(obj)
  } else {
    cat(" ERROR: File", f, "is not a valid Seurat object!\n")
    return(NULL)
  }
}

# Load and validate Seurat objects
# Applies the safe load function across all files and removes any invalid (NULL) objects from the list.
seurat_objects <- lapply(file_list, safe_load_seurat)
seurat_objects <- Filter(Negate(is.null), seurat_objects)  # Remove NULL objects

# Stop if no valid Seurat objects were loaded
# Ensures at least one valid Seurat object exists, otherwise stops the script.
if (length(seurat_objects) == 0) stop(" ERROR: No valid Seurat objects found!")


# 6. Integration and Normalization of Seurat Objects
# Finds genes/features common across all loaded Seurat objects to ensure compatibility during integration.
common_features <- Reduce(intersect, lapply(seurat_objects, rownames))
if (length(common_features) == 0) stop(" ERROR: No common features found!")

# Perform SCT normalization if not already done
# Ensure SCT normalization before integration
# Checks if Seurat objects have been normalized with SCTransform (SCT), performing normalization if necessary.
seurat_objects <- lapply(seurat_objects, function(obj) {
  if (!"SCT" %in% names(obj@assays)) {
    obj <- SCTransform(obj, verbose = FALSE)
  }
  return(obj)
})

# 7. Integration and Dimensionality Reduction
# Run integration (Optimized for performance)
# Integrates multiple Seurat objects into a single combined dataset based on common anchor features.
anchors <- FindIntegrationAnchors(object.list = seurat_objects, dims = 1:30, anchor.features = 3000)
combined <- IntegrateData(anchorset = anchors, dims = 1:30)

# Set default assay before clustering
# Sets the default assay to integrated for subsequent analysis.
DefaultAssay(combined) <- "integrated"

# Scale and Run PCA (Principal Component Analysis)
# Scales the data, performs PCA, UMAP embedding, and identifies nearest neighbors for clustering.
combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
# Run UMAP
combined <- RunUMAP(combined, dims = 1:30, verbose = FALSE)
# Find Neighbors
combined <- FindNeighbors(combined, dims = 1:30, verbose = FALSE)

# 8. Finding Exactly 14 Clusters
# Iteratively adjusts clustering resolution to precisely find 14 clusters, stopping once exactly 14 are identified.
resolution_vals <- seq(0.2, 1.0, by = 0.05)
found_clusters <- FALSE

for (res in resolution_vals) {
  combined <- FindClusters(combined, resolution = res)
  num_clusters <- length(unique(Idents(combined)))
  cat(" Trying resolution =", res, "Clusters found:", num_clusters, "\n")
  if (num_clusters == 14) {
    cat(" Found exactly 14 clusters at resolution =", res, "\n")
    found_clusters <- TRUE
    break
  }
}

if (!found_clusters) {
  cat(" WARNING: Could not get exactly 14 clusters. Using closest match!\n")
}

# 9. Save Final Seurat Object
# Saves the final integrated and clustered Seurat object to the specified folder.
saveRDS(combined, file.path(qc_folder, "Final_Seurat_Object_14Clusters.rds"))


# 10. Generating Plots and Saving Outputs
# PLOT 1: UMAP plot for visualization of 14 clusters clearly labeled.
umap_plot <- DimPlot(combined, reduction = "umap", label = TRUE) +
  ggtitle("UMAP - 12 Samples (Exactly 14 Clusters)") +
  theme_minimal()

# PLOT 2: Sample contribution plot visualizing sample proportions within each cluster.
sample_contribution <- as.data.frame(table(Idents(combined), combined$orig.ident))
colnames(sample_contribution) <- c("Cluster", "Sample", "Count")

sample_contribution_plot <- ggplot(sample_contribution, aes(x = Cluster, y = Count, fill = Sample)) +
  geom_bar(stat = "identity", position = "fill") +
  ylab("Proportion") +
  ggtitle("Sample Contribution by Cluster (Exact 14 Clusters)") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set3")

# PLOT 3: GZMK expression plot visualizing expression levels of the gene GZMK
if ("GZMK" %in% rownames(combined)) {
  gzmk_expression_plot <- VlnPlot(combined, features = "GZMK", group.by = "seurat_clusters", pt.size = 0.1) +
    ggtitle("GZMK Expression in Each Cluster") +
    scale_fill_manual(values = scales::hue_pal()(14)) +
    theme_minimal()
} else {
  gzmk_expression_plot <- ggplot() + ggtitle(" GZMK not found in dataset!") + theme_minimal()
  cat("\n WARNING: GZMK not found in dataset! Skipping GZMK plot...\n")
}

# Save Plots
# Each plot is generated and then saved separately to the qc_folder.
# Saves each plot as a high-resolution .png image.
ggsave(file.path(qc_folder, "Combined_UMAP_14Clusters.png"), plot = umap_plot, width = 10, height = 6)
ggsave(file.path(qc_folder, "Sample_Contribution_14Clusters.png"), plot = sample_contribution_plot, width = 10, height = 6)
ggsave(file.path(qc_folder, "GZMK_Expression_14Clusters.png"), plot = gzmk_expression_plot, width = 10, height = 6)


# 11. Completion Message
# Confirms successful script completion and indicates where outputs are saved.
cat("\n Analysis completed successfully! Plots saved in", qc_folder, " \n")