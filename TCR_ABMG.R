# Shamim_Roshni - Full Pipeline with Additional Plots

# Setup ---
options(future.globals.maxSize = 450 * 1024^3)
options(repos = c(CRAN = "https://cran.r-project.org"))

# Add this to make sure all packages are found:
.libPaths(c("/home/abmshamimuh/R/x86_64-pc-linux-gnu-library/4.4", .libPaths()))


# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(readxl)
  library(stringr)
})


# Define paths
data_folder <- "/home/abmshamimuh/"
qc_folder <- file.path(data_folder, "QC_Results4")
tcr_file <- file.path(data_folder, "Corridoni TCRs.xlsx")
seurat_rds <- file.path(qc_folder, "Final_Seurat_Object_14Clusters.rds")

# Load Seurat object
combined <- readRDS(seurat_rds)

# Load TCR metadata
tcr_data <- read_excel(tcr_file, sheet = "counts with matches")
tcr_data <- tcr_data %>% rename_with(make.names)  # handle duplicate column names
tcr_data$barcode_clean <- str_replace(tcr_data$barcode, "-1$", "")
tcr_data$viral_status <- ifelse(is.na(tcr_data$organism), "Non-Viral", "Viral")

# Extract first 16 bases from Seurat cell names as barcodes
combined$seurat_barcode <- substr(Cells(combined), 1, 16)

# Debug output to check
cat("Example matched Seurat barcode:\n")
print(head(combined$seurat_barcode, 5))
cat("\n Example TCR barcode:\n")
print(head(tcr_data$barcode_clean, 5))

# Join metadata
barcode_df <- data.frame(seurat_barcode = combined$seurat_barcode, row.names = Cells(combined))
barcode_df <- left_join(barcode_df, tcr_data[, c("barcode_clean", "viral_status")],
                        by = c("seurat_barcode" = "barcode_clean"))
barcode_df$viral_status[is.na(barcode_df$viral_status)] <- "Non-Viral"
combined$Viral <- barcode_df$viral_status

# Confirm match
cat("\n Viral TCRs matched in Seurat object:", sum(combined$Viral == "Viral"), "\n")





# 1. Split viral vs non-viral
viral_cells <- subset(combined, subset = Viral == "Viral")
nonviral_cells <- subset(combined, subset = Viral == "Non-Viral")

# 2. Extract UMAP and add a "Group" label
df_nonviral <- Embeddings(nonviral_cells, reduction = "umap") %>%
  as.data.frame() %>%
  rename(UMAP_1 = umap_1, UMAP_2 = umap_2) %>%
  mutate(Group = "Non-Viral TCRs")

df_viral <- Embeddings(viral_cells, reduction = "umap") %>%
  as.data.frame() %>%
  rename(UMAP_1 = umap_1, UMAP_2 = umap_2) %>%
  mutate(Group = "Viral TCRs")

# 3. Combine data
df_all <- rbind(df_nonviral, df_viral)

# 4. Plot with custom legend


# Count number of viral-labeled cells in UMAP
n_viral <- sum(df_all$Group == "Viral TCRs")

# 4. Plot with custom legend and viral count in title
viral_umap <- ggplot(df_all, aes(x = UMAP_1, y = UMAP_2, color = Group)) +
  geom_point(size = 0.4) +
  scale_color_manual(values = c("Non-Viral TCRs" = "gray80", "Viral TCRs" = "red")) +
  labs(
    title = paste0("Combined UMAP (12 Samples) with Viral TCR Overlay (n = ", n_viral, " Viral Cells)\nForced 14 Clusters"),
    x = "UMAP 1", y = "UMAP 2", color = "TCR Type"
  ) +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.6),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

# 5. Save final plot
ggsave(file.path(qc_folder, "Plot4_UMAP_Viral_vs_NonViral_Overlay_CIRCLES_LABELED_AXIS_LEGEND.png"),
       plot = viral_umap, width = 10, height = 7, dpi = 300)

# 6. View it
print(viral_umap)


# For 371
# Load necessary libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(readxl)
library(stringr)

# Set paths
qc_folder <- "/home/abmshamimuh/QC_Results4"
tcr_file <- "/home/abmshamimuh/Corridoni TCRs.xlsx"
seurat_file <- file.path(qc_folder, "Final_Seurat_Object_14Clusters.rds")

# Load Seurat object
combined <- readRDS(seurat_file)

# Load TCR metadata and clean
tcr_data <- read_excel(tcr_file, sheet = "counts with matches") %>%
  rename_with(make.names)
tcr_data$barcode_clean <- str_replace(tcr_data$barcode, "-1$", "")
tcr_data$viral_status <- ifelse(is.na(tcr_data$organism), "Non-Viral", "Viral")

# Clean Seurat barcodes and join with metadata
combined$seurat_barcode <- substr(Cells(combined), 1, 16)
barcode_df <- data.frame(seurat_barcode = combined$seurat_barcode, row.names = Cells(combined))
barcode_df <- left_join(barcode_df, tcr_data[, c("barcode_clean", "viral_status", "cdr3_beta")],
                        by = c("seurat_barcode" = "barcode_clean"))
barcode_df$viral_status[is.na(barcode_df$viral_status)] <- "Non-Viral"
combined$Viral <- barcode_df$viral_status
combined$CDR3beta <- barcode_df$cdr3_beta

# Identify all viral-labeled Seurat cells
viral_cells <- subset(combined, subset = Viral == "Viral")

# Create a DataFrame with barcode and cell name
viral_info <- data.frame(
  barcode = viral_cells$seurat_barcode,
  cell = Cells(viral_cells),
  row.names = NULL
)

# Deduplicate by barcode: take the first cell for each unique barcode
unique_viral_cells <- viral_info %>%
  distinct(barcode, .keep_all = TRUE)  # 1 row per unique barcode

# Double check count
n_viral_unique <- nrow(unique_viral_cells)  # Should be 378

# UMAP DataFrame with labels
df_all <- Embeddings(combined, "umap") %>%
  as.data.frame() %>%
  mutate(
    UMAP_1 = .[[1]],
    UMAP_2 = .[[2]],
    Viral_Label = ifelse(rownames(.) %in% unique_viral_cells$cell, "Viral TCRs", "Non-Viral TCRs")
  )

table(df_all$Viral_Label)

# Count for title
n_viral_unique <- sum(df_all$Viral_Label == "Viral TCRs")

# Plot
viral_umap <- ggplot(df_all, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(data = subset(df_all, Viral_Label == "Non-Viral TCRs"),
             aes(color = Viral_Label), size = 0.3, alpha = 0.5) +
  geom_point(data = subset(df_all, Viral_Label == "Viral TCRs"),
             aes(color = Viral_Label), size = 0.9, alpha = 0.9) +
  scale_color_manual(values = c("Non-Viral TCRs" = "gray80", "Viral TCRs" = "red")) +
  labs(
    title = paste0("UMAP: Unique Viral TCRs in Red (n = ", n_viral_unique, ")\nOne Cell per Barcode"),
    x = "UMAP 1", y = "UMAP 2", color = "TCR Type"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10)
  )

# Save
ggsave(
  file.path(qc_folder, "371Plot1_UMAP_Unique_Viral_TCRs_Final.png"),
  plot = viral_umap, width = 10, height = 7, dpi = 300
)

# Print to screen
print(viral_umap)


# Plot i) Unique Viral TCRs by Cluster (Percentage)
cluster_data <- data.frame(
  Cluster = combined$seurat_clusters[unique_viral_cells$cell]
) %>%
  count(Cluster, name = "Count") %>%
  mutate(Percentage = round(100 * Count / sum(Count), 1))

plot_cluster <- ggplot(cluster_data, aes(x = Cluster, y = Percentage, fill = Cluster)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(Percentage, "%")), vjust = -0.3, size = 4.5) +
  labs(
    title = paste0("Unique Viral TCRs by Cluster (n = ", nrow(unique_viral_cells), ")"),
    x = "Cluster", y = "Percentage of Unique Viral TCRs"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none"
  )

ggsave(file.path(qc_folder, "371Plot2_Unique_ViralTCR_Percent_by_Cluster.png"),
       plot = plot_cluster, width = 8, height = 6, dpi = 300)


# Plot ii) Unique Viral TCRs by Sample (Percentage)
sample_data <- data.frame(
  Sample = combined$orig.ident[unique_viral_cells$cell]
) %>%
  count(Sample, name = "Count") %>%
  mutate(Percentage = round(100 * Count / sum(Count), 1))

plot_sample <- ggplot(sample_data, aes(x = Sample, y = Percentage, fill = Sample)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(Percentage, "%")), vjust = -0.3, size = 4.5) +
  labs(
    title = paste0("Unique Viral TCRs by Sample (n = ", nrow(unique_viral_cells), ")"),
    x = "Sample", y = "Percentage of Unique Viral TCRs"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(file.path(qc_folder, "371Plot3_Unique_ViralTCR_Percent_by_Sample.png"),
       plot = plot_sample, width = 9, height = 6, dpi = 300)


# Add Seurat cluster info to unique viral cells
unique_viral_cells$cluster <- combined$seurat_clusters[unique_viral_cells$cell]

# Add GZMK expression info
gzmk_expr <- FetchData(combined, vars = "GZMK")[unique_viral_cells$cell, , drop = FALSE]
unique_viral_cells$GZMK <- gzmk_expr$GZMK
unique_viral_cells$GZMK_Pos <- ifelse(unique_viral_cells$GZMK > 0, 1, 0)

# Count % of GZMK+ cells per cluster
gzmk_by_cluster <- unique_viral_cells %>%
  group_by(cluster) %>%
  summarise(
    Total = n(),
    GZMK_Pos = sum(GZMK_Pos),
    Percent = round(GZMK_Pos / Total * 100, 1)
  )

# Plot
p_cluster <- ggplot(gzmk_by_cluster, aes(x = factor(cluster), y = Percent, fill = factor(cluster))) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(Percent, "%")), vjust = -0.3, size = 4) +
  labs(
    title = paste0("GZMK+ Cells in Unique Viral TCRs by Cluster (n = ", nrow(unique_viral_cells), ")"),
    x = "Seurat Cluster", y = "GZMK+ Percentage"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "none"
  )

ggsave(file.path(qc_folder, "371Plot4_GZMK_UniqueViral_by_Cluster_Percent.png"), plot = p_cluster, width = 8, height = 6)


# statistical tests

# Libraries (already loaded in your code context)
library(dplyr)
library(Seurat)
library(ggplot2)

# Step 1: Add cluster and GZMK info to unique viral TCR cells
unique_viral_cells$cluster <- combined$seurat_clusters[unique_viral_cells$cell]
gzmk_expr <- FetchData(combined, vars = "GZMK")[unique_viral_cells$cell, , drop = FALSE]
unique_viral_cells$GZMK <- gzmk_expr$GZMK
unique_viral_cells$GZMK_Pos <- ifelse(unique_viral_cells$GZMK > 0, 1, 0)

# Step 2: Summarize % of GZMK+ per cluster
gzmk_by_cluster <- unique_viral_cells %>%
  group_by(cluster) %>%
  summarise(
    Total = n(),
    GZMK_Pos = sum(GZMK_Pos),
    Percent = round(GZMK_Pos / Total * 100, 1)
  )

# Step 3: Create contingency table
fisher_table <- table(unique_viral_cells$cluster, unique_viral_cells$GZMK_Pos)

# ---- Chi-Square test across all clusters ----
chisq_result <- chisq.test(fisher_table)
cat("371Plot4_Chi-Square Test (All Clusters): p =", format.pval(chisq_result$p.value, digits = 10, scientific = TRUE), "\n")

# ---- Fisher’s Exact test (compare 2 clusters only, e.g., 10 vs 3) ----
fisher_result <- fisher.test(fisher_table[c("4", "6"), ])
cat("371Plot4_Fisher’s Exact Test (Cluster 4 vs 6): p =", signif(fisher_result$p.value, 4), "\n")


# save the results

# Step 3: Create contingency table
fisher_table <- table(unique_viral_cells$cluster, unique_viral_cells$GZMK_Pos)

# ---- Chi-Square test across all clusters (0–13) ----
chisq_result <- chisq.test(fisher_table)
chisq_pval <- format.pval(chisq_result$p.value, digits = 10, scientific = TRUE)

# ---- Fisher’s Exact Test between Cluster 4 and Cluster 6 ----
fisher_result <- fisher.test(fisher_table[c("4", "6"), ])
fisher_pval <- signif(fisher_result$p.value, 4)

# ---- Prepare Results Text ----
results_text <- paste0(
  "Statistical Test Results for GZMK+ Distribution by Cluster (371 Unique Viral TCR Cells)\n",
  "----------------------------------------------------------------------------------------\n",
  "Chi-Square Test (All 14 Clusters):\n",
  "  p-value = ", chisq_pval, "\n\n",
  "Fisher’s Exact Test (Cluster 4 vs Cluster 6):\n",
  "  p-value = ", fisher_pval, "\n"
)

# ---- Save to file ----
writeLines(results_text, "/home/abmshamimuh/QC_Results4/371Plot4_StatisticalTests.txt")

# ---- Also print to console ----
cat(results_text)






# Add sample info to unique viral cells
unique_viral_cells$sample <- combined$orig.ident[unique_viral_cells$cell]

# Count % of GZMK+ cells per sample
gzmk_by_sample <- unique_viral_cells %>%
  group_by(sample) %>%
  summarise(
    Total = n(),
    GZMK_Pos = sum(GZMK_Pos),
    Percent = round(GZMK_Pos / Total * 100, 1)
  )

# Plot
p_sample <- ggplot(gzmk_by_sample, aes(x = sample, y = Percent, fill = sample)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(Percent, "%")), vjust = -0.3, size = 4) +
  labs(
    title = paste0("GZMK+ Cells in Unique Viral TCRs by Sample (n = ", nrow(unique_viral_cells), ")"),
    x = "Sample", y = "GZMK+ Percentage"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "none"
  )

ggsave(file.path(qc_folder, "371Plot5_GZMK_UniqueViral_by_Sample_Percent.png"), plot = p_sample, width = 10, height = 6)

# Subset Seurat object to only 371 unique viral TCRs (actual present barcodes)
viral_seurat_cells <- unique_viral_cells$cell  # 371 cells

# UMAP with GZMK expression
gzm_expr <- FetchData(combined, vars = c("GZMK"), cells = viral_seurat_cells)
umap_coords <- Embeddings(combined, "umap")[viral_seurat_cells, ]
df_gzmk_umap <- data.frame(UMAP_1 = umap_coords[, 1],
                           UMAP_2 = umap_coords[, 2],
                           GZMK = gzm_expr$GZMK)

# Plot UMAP expression
gzmk_umap <- ggplot(df_gzmk_umap, aes(x = UMAP_1, y = UMAP_2, color = GZMK)) +
  geom_point(size = 1.2) +
  scale_color_viridis_c(option = "plasma") +
  labs(title = paste0("UMAP: GZMK Expression in Unique Viral TCRs (n = ", nrow(df_gzmk_umap), ")"),
       x = "UMAP 1", y = "UMAP 2", color = "GZMK") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Save UMAP
ggsave(file.path(qc_folder, "371Plot6_GZMK_UMAP_UniqueViralTCRs.png"), gzmk_umap, width = 10, height = 7, dpi = 300)
print(gzmk_umap)

# Label all cells as Viral or Non-Viral
viral_flags <- ifelse(Cells(combined) %in% unique_viral_cells$cell, "Viral", "Non-Viral")
combined$Viral_GZMK_Group <- factor(viral_flags, levels = c("Non-Viral", "Viral"))

# Violin plot of GZMK
vln_plot <- VlnPlot(combined,
                    features = "GZMK",
                    group.by = "Viral_GZMK_Group",
                    pt.size = 0.1) +
  labs(title = paste0("GZMK Expression in Viral vs Non-Viral Cells (n = ", nrow(unique_viral_cells), " Viral)")) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Save
ggsave(file.path(qc_folder, "371Plot7_Violin_GZMK_Viral_vs_NonViral.png"), vln_plot, width = 6, height = 5, dpi = 300)
print(vln_plot)

# 1. Extract the 371 unique Seurat cell names matched to viral TCRs
viral_seurat_cells <- unique_viral_cells$cell  # Already deduplicated

# 2. Add Viral_Label column to all Seurat cells
combined$Viral_Label <- ifelse(Cells(combined) %in% viral_seurat_cells, "Viral", "Non-Viral")

# 3. Plot violin plot of GZMK (Viral vs Non-Viral — using deduplicated labels)
vln_plot <- VlnPlot(combined, features = "GZMK", group.by = "Viral_Label", pt.size = 0.1) +
  labs(title = paste0("GZMK Expression in Viral vs Non-Viral Cells\n(Using 1 Cell per Barcode, n = ", length(viral_seurat_cells), " Viral)")) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Save
ggsave(file.path(qc_folder, "371Plot8_Violin_GZMK_371Unique_Viral_vs_NonViral.png"),
       plot = vln_plot, width = 6, height = 5, dpi = 300)

# Print
print(vln_plot)



# Load required libraries
library(dplyr)
library(ggplot2)
library(readxl)
library(stringr)
library(RColorBrewer)

# File paths
tcr_file <- "/home/abmshamimuh/Corridoni TCRs.xlsx"
output_dir <- "/home/abmshamimuh/QC_Results4"

# Load and clean TCR metadata
tcr_data <- read_excel(tcr_file, sheet = "counts with matches") %>%
  rename_with(make.names)

tcr_data$barcode_clean <- str_replace(tcr_data$barcode, "-1$", "")
tcr_data$viral_status <- ifelse(is.na(tcr_data$organism), "Non-Viral", "Viral")

# Assign group (HC or UC)
tcr_data$group <- ifelse(tcr_data$Sample %in% c("S21", "S22", "S23"), "HC",
                         ifelse(tcr_data$Sample %in% c("S24", "S33", "S34"), "UC", NA))

# Filter to valid viral entries only
viral_tcr_data <- tcr_data %>%
  filter(viral_status == "Viral", !is.na(group), !is.na(organism), organism != "") %>%
  mutate(organism = trimws(organism)) %>%
  filter(tolower(organism) != "na")

# Set color palette for organisms
all_organisms <- sort(unique(viral_tcr_data$organism))
viral_tcr_data$organism <- factor(viral_tcr_data$organism, levels = all_organisms)
set.seed(42)
organism_colors <- setNames(colorRampPalette(brewer.pal(8, "Set3"))(length(all_organisms)), all_organisms)

# Pie Chart for HC
hc_data <- viral_tcr_data %>% filter(group == "HC")
if (nrow(hc_data) > 0) {
  plot_7 <- hc_data %>%
    count(organism) %>%
    mutate(perc = 100 * n / sum(n),
           label = ifelse(perc > 1, paste0(round(perc, 1), "%"), "")) %>%
    ggplot(aes(x = "", y = perc, fill = organism)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar("y") +
    geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 4) +
    scale_fill_manual(values = organism_colors) +
    labs(title = paste0("Viral Specific TCRs by HC (n = ", sum(hc_data$organism != ""), ")")) +
    theme_void() +
    theme(plot.title = element_text(face = "bold"))

  ggsave(file.path(output_dir, "280Plot7_Viral_HC_Organism_Perc.png"), plot_7, width = 7, height = 6)
}

# Pie Chart for UC
uc_data <- viral_tcr_data %>% filter(group == "UC")
if (nrow(uc_data) > 0) {
  plot_8 <- uc_data %>%
    count(organism) %>%
    mutate(perc = 100 * n / sum(n),
           label = ifelse(perc > 1, paste0(round(perc, 1), "%"), "")) %>%
    ggplot(aes(x = "", y = perc, fill = organism)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar("y") +
    geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 4) +
    scale_fill_manual(values = organism_colors) +
    labs(title = paste0("Viral Specific TCRs by UC (n = ", sum(uc_data$organism != ""), ")")) +
    theme_void() +
    theme(plot.title = element_text(face = "bold"))

  ggsave(file.path(output_dir, "280Plot8_Viral_UC_Organism_Perc.png"), plot_8, width = 7, height = 6)
}



# Load required libraries
library(Seurat)
library(dplyr)
library(readxl)
library(stringr)

# File paths
tcr_file <- "/home/abmshamimuh/Corridoni TCRs.xlsx"
seurat_file <- "/home/abmshamimuh/QC_Results4/Final_Seurat_Object_14Clusters.rds"

# Load Seurat object
combined <- readRDS(seurat_file)
combined$seurat_barcode <- substr(Cells(combined), 1, 16)

# Load and clean TCR metadata
tcr_data <- read_excel(tcr_file, sheet = "counts with matches") %>%
  rename_with(make.names)

tcr_data$barcode_clean <- str_replace(tcr_data$barcode, "-1$", "")
tcr_data$viral_status <- ifelse(is.na(tcr_data$organism), "Non-Viral", "Viral")

#  Annotate group for HC vs UC (6 labeled samples)
tcr_data$group <- ifelse(tcr_data$Sample %in% c("S21", "S22", "S23"), "HC",
                         ifelse(tcr_data$Sample %in% c("S24", "S33", "S34"), "UC", NA))

#  Step 1: Filter 280 unique viral TCR barcodes from 6 matched samples
tcr_viral_280 <- tcr_data %>%
  filter(group %in% c("HC", "UC"),
         viral_status == "Viral",
         !is.na(organism),
         organism != "",
         tolower(organism) != "na") %>%
  distinct(barcode_clean, .keep_all = TRUE)

barcodes_280 <- tcr_viral_280$barcode_clean
cat(" Number of unique viral barcodes in 280 subset:", length(barcodes_280), "\n")

# Step 2: Assign Viral status into Seurat object
barcode_df <- data.frame(seurat_barcode = combined$seurat_barcode, row.names = Cells(combined))
barcode_df <- left_join(barcode_df, tcr_data[, c("barcode_clean", "viral_status")],
                        by = c("seurat_barcode" = "barcode_clean"))
barcode_df$viral_status[is.na(barcode_df$viral_status)] <- "Non-Viral"
combined$Viral <- barcode_df$viral_status

#  Step 3: Get 371 barcodes from Seurat object
seurat_barcodes_371 <- unique(combined$seurat_barcode[combined$Viral == "Viral"])
cat(" Total unique viral barcodes in Seurat (371 set):", length(seurat_barcodes_371), "\n")

# Step 4: Compare overlap
shared_barcodes <- intersect(barcodes_280, seurat_barcodes_371)
cat("Overlap between 280 and 371 sets:", length(shared_barcodes), "\n")

# Differences
only_in_280 <- setdiff(barcodes_280, seurat_barcodes_371)
only_in_371 <- setdiff(seurat_barcodes_371, barcodes_280)
cat(" Barcodes in 280 not in Seurat:", length(only_in_280), "\n")
cat(" Barcodes in Seurat not in 280:", length(only_in_371), "\n")

# Optional: save unmatched
writeLines(only_in_280, "/home/abmshamimuh/QC_Results4/barcodes_280_not_in_seurat.txt")
writeLines(only_in_371, "/home/abmshamimuh/QC_Results4/barcodes_371_not_in_280.txt")

# Step 1: Count unique viral TCRs in the 280 subset (Excel metadata)
barcodes_280 <- tcr_data %>%
  filter(Sample %in% c("S21", "S22", "S23", "S24", "S33", "S34"),
         viral_status == "Viral",
         !is.na(organism),
         organism != "",
         tolower(organism) != "na") %>%
  pull(barcode_clean) %>%
  unique()

cat(" Unique viral TCR barcodes in 280 sample subset:", length(barcodes_280), "\n")

# Step 2: Get viral barcodes from Seurat
seurat_barcodes_371 <- combined$seurat_barcode[combined$Viral == "Viral"] %>% unique()

# Step 3: Intersect
barcodes_280_in_seurat <- intersect(barcodes_280, seurat_barcodes_371)

cat(" Barcodes from 280 that are present in Seurat (part of 371):", length(barcodes_280_in_seurat), "\n")

writeLines(barcodes_280_in_seurat, file.path(qc_folder, "Barcodes_280_matched_to_Seurat.txt"))
writeLines(setdiff(barcodes_280, seurat_barcodes_371), file.path(qc_folder, "Barcodes_280_NOT_in_Seurat.txt"))
writeLines(setdiff(seurat_barcodes_371, barcodes_280), file.path(qc_folder, "Barcodes_371_NOT_in_280.txt"))


# 280 are also Unique viral TCRs, 30 are too

# "Violin and bar plots showing GZMK expression among 30 unique viral TCR-bearing Seurat cells (from matched HC and UC samples). 
# Only barcodes with exact cell matches were used to avoid clonal bias.
# Exactly 30 unique Seurat cells based on 280 HC/UC viral barcodes (1 cell per barcode, no clones).
# Only 30 of these 280 barcodes actually exist in the Seurat object with 371
# If you're analyzing GZMK expression by HC vs UC, only 30 Seurat cells (from the 280 set) are actually usable.


# Load libraries
library(Seurat)
library(dplyr)
library(readxl)
library(stringr)
library(ggplot2)

# File paths
tcr_file <- "/home/abmshamimuh/Corridoni TCRs.xlsx"
seurat_file <- "/home/abmshamimuh/QC_Results4/Final_Seurat_Object_14Clusters.rds"

# Load Seurat object
combined <- readRDS(seurat_file)
combined$barcode_clean <- substr(Cells(combined), 1, 16)

# Load and clean TCR metadata
tcr_data <- read_excel(tcr_file, sheet = "counts with matches") %>%
  rename_with(make.names)
tcr_data$barcode_clean <- str_replace(tcr_data$barcode, "-1$", "")
tcr_data$viral_status <- ifelse(is.na(tcr_data$organism), "Non-Viral", "Viral")
tcr_data$group <- case_when(
  tcr_data$Sample %in% c("S21", "S22", "S23") ~ "HC",
  tcr_data$Sample %in% c("S24", "S33", "S34") ~ "UC",
  TRUE ~ NA_character_
)

# Step 1: Filter clean 280 and find overlap with Seurat
tcr_280 <- tcr_data %>%
  filter(group %in% c("HC", "UC"),
         viral_status == "Viral",
         !is.na(organism),
         organism != "",
         tolower(organism) != "na") %>%
  distinct(barcode_clean, .keep_all = TRUE)

barcodes_280 <- tcr_280$barcode_clean
seurat_barcodes <- unique(substr(Cells(combined), 1, 16))
matched_30 <- intersect(barcodes_280, seurat_barcodes)

# Step 2: Match one cell per barcode
cell_df <- data.frame(cell = Cells(combined), barcode_clean = substr(Cells(combined), 1, 16))
selected_cells <- cell_df %>%
  filter(barcode_clean %in% matched_30) %>%
  group_by(barcode_clean) %>%
  slice(1) %>%
  ungroup()

# Step 3: Annotate Seurat object
group_map <- tcr_280 %>%
  filter(barcode_clean %in% matched_30) %>%
  select(barcode_clean, group) %>%
  distinct()

group_meta <- left_join(selected_cells, group_map, by = "barcode_clean")
combined$Group30_exact <- NA
combined$Group30_exact[group_meta$cell] <- group_meta$group

# Step 4: Subset Seurat to 30 cells
cells_30_exact <- group_meta$cell
combined_30 <- subset(combined, cells = cells_30_exact)

# ---- Plot 1: GZMK Expression Violin ----
p1 <- VlnPlot(combined_30,
              features = "GZMK",
              group.by = "Group30_exact",
              pt.size = 0.5) +
  labs(title = "GZMK Expression in HC vs UC (n = 30 Unique Viral TCRs)",
       x = "Group", y = "GZMK Expression Level") +
  scale_fill_manual(values = c("HC" = "salmon", "UC" = "turquoise3")) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    axis.title = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  )

ggsave("/home/abmshamimuh/QC_Results4/Plot13_GZMK_HC_UC_30cells_Manuscript.png", p1, width = 7.5, height = 6)

# ---- Plot 2: GZMK+ Percentage Barplot ----
expr_data <- FetchData(combined_30, vars = c("GZMK", "Group30_exact")) %>%
  mutate(GZMK_Positive = GZMK > 0)

gzmk_summary <- expr_data %>%
  group_by(Group30_exact) %>%
  summarise(
    Total_Cells = n(),
    GZMK_Positive = sum(GZMK_Positive),
    Percent_Positive = round(100 * GZMK_Positive / Total_Cells, 1)
  )

p2 <- ggplot(gzmk_summary, aes(x = Group30_exact, y = Percent_Positive, fill = Group30_exact)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = paste0(Percent_Positive, "%")), vjust = -0.5, size = 5) +
  scale_fill_manual(values = c("HC" = "salmon", "UC" = "turquoise3")) +
  ylim(0, 100) +
  labs(
    title = "Percentage of GZMK+ Cells in HC vs UC (n = 30 Unique Viral TCRs)",
    x = "Group", y = "GZMK+ Cells (%)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    axis.title = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 12),
    legend.position = "none"
  )

ggsave("/home/abmshamimuh/QC_Results4/Plot14_GZMK_Percent_HC_UC_30cells_Manuscript.png", p2, width = 6.5, height = 5)

# Print to scree
print(p1)
print(p2)
print(gzmk_summary)


# Fetch GZMK expression and HC/UC labels
expr_data <- FetchData(combined_30, vars = c("GZMK", "Group30_exact")) %>%
  mutate(GZMK_Positive = GZMK > 0)

# ------------------------
# Wilcoxon Rank-Sum Test
# ------------------------
wilcox_result <- wilcox.test(GZMK ~ Group30_exact, data = expr_data)
cat("Wilcoxon Rank-Sum Test p-value:", wilcox_result$p.value, "\n")

# ------------------------
# Fisher’s Exact Test
# ------------------------
# Contingency table: group × GZMK_Positive
fisher_table <- table(expr_data$Group30_exact, expr_data$GZMK_Positive)
fisher_result <- fisher.test(fisher_table)
cat("Fisher’s Exact Test p-value:", fisher_result$p.value, "\n")


# save results
# ------------------------
# Wilcoxon Rank-Sum Test
# ------------------------
wilcox_result <- wilcox.test(GZMK ~ Group30_exact, data = expr_data)

# ------------------------
# Fisher’s Exact Test
# ------------------------
fisher_table <- table(expr_data$Group30_exact, expr_data$GZMK_Positive)
fisher_result <- fisher.test(fisher_table)

# ------------------------
# Save results to file
# ------------------------
results_text <- paste0(
  "Statistical Test Results for GZMK Expression (Plot13–14)\n",
  "--------------------------------------------------------\n",
  "Wilcoxon Rank-Sum Test:\n",
  "  p-value = ", signif(wilcox_result$p.value, 4), "\n\n",
  "Fisher’s Exact Test:\n",
  "  p-value = ", signif(fisher_result$p.value, 4), "\n"
)

# Write to file
writeLines(results_text, "/home/abmshamimuh/QC_Results4/Plot13_14_StatisticalTests.txt")

# Also print to console
cat(results_text)





# ------------------------
# Violin Plot with Annotated p-values
# ------------------------
p1 <- VlnPlot(combined_30,
              features = "GZMK",
              group.by = "Group30_exact",
              pt.size = 0.5) +
  labs(title = "GZMK Expression in HC vs UC (n = 30 Unique Viral TCRs)",
       x = "Group", y = "GZMK Expression Level") +
  scale_fill_manual(values = c("HC" = "salmon", "UC" = "turquoise3")) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    axis.title = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  ) +
  annotate("text",
           x = 1.5,
           y = max(expr_data$GZMK, na.rm = TRUE) * 1.1,
           label = paste0("Wilcoxon p = ", signif(wilcox_result$p.value, 3)),
           size = 5,
           fontface = "italic") +
  annotate("text",
           x = 1.5,
           y = max(expr_data$GZMK, na.rm = TRUE) * 1.2,
           label = paste0("Fisher's p = ", signif(fisher_result$p.value, 3)),
           size = 5,
           fontface = "italic")

# Save and print
ggsave("/home/abmshamimuh/QC_Results4/p-valuePlot13_GZMK_HC_UC_30cells_Annotated.png", p1, width = 7.5, height = 6)
print(p1)










# From 378 cell/ but it will 371 cell in plot, because we have 371 unique in seruet, 7 missing 

# Load required libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(readxl)
library(stringr)

# Paths
qc_folder <- "/home/abmshamimuh/QC_Results4"
tcr_file <- "/home/abmshamimuh/Corridoni TCRs.xlsx"
seurat_file <- file.path(qc_folder, "Final_Seurat_Object_14Clusters.rds")

# Load Seurat object
combined <- readRDS(seurat_file)

# Load and clean TCR metadata
tcr_data <- read_excel(tcr_file, sheet = "counts with matches") %>%
  rename_with(make.names)
tcr_data$barcode_clean <- str_replace(tcr_data$barcode, "-1$", "")
tcr_data$viral_status <- ifelse(is.na(tcr_data$organism), "Non-Viral", "Viral")

# Extract 378 unique viral barcodes from Excel
excel_viral_barcodes <- unique(tcr_data$barcode_clean[tcr_data$viral_status == "Viral"])
length(excel_viral_barcodes)  # Should be 378

# Match Seurat cells to Excel viral barcodes
combined$seurat_barcode <- substr(Cells(combined), 1, 16)
viral_cells <- subset(combined, subset = seurat_barcode %in% excel_viral_barcodes)

# Deduplicate: Take 1 cell per Excel barcode (even if barcode repeated)
viral_df <- data.frame(
  barcode = viral_cells$seurat_barcode,
  cell = Cells(viral_cells),
  stringsAsFactors = FALSE
)
viral_unique_df <- viral_df %>%
  distinct(barcode, .keep_all = TRUE)  # 1 row per barcode from Excel

# Confirm total viral cells selected for UMAP
cat(" Viral cells used for UMAP overlay:", nrow(viral_unique_df), "\n")  # Expect 378

# Build UMAP overlay labels
df_all <- Embeddings(combined, "umap") %>%
  as.data.frame() %>%
  mutate(
    UMAP_1 = .[[1]],
    UMAP_2 = .[[2]],
    Viral_Label = ifelse(rownames(.) %in% viral_unique_df$cell, "Viral TCRs", "Non-Viral TCRs")
  )

# Plot UMAP
umap_plot <- ggplot(df_all, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(data = subset(df_all, Viral_Label == "Non-Viral TCRs"),
             aes(color = Viral_Label), size = 0.3, alpha = 0.5) +
  geom_point(data = subset(df_all, Viral_Label == "Viral TCRs"),
             aes(color = Viral_Label), size = 0.9, alpha = 0.9) +
  scale_color_manual(values = c("Non-Viral TCRs" = "gray80", "Viral TCRs" = "red")) +
  labs(
    title = paste0("UMAP Overlay: 378 Viral TCRs (One per Excel Barcode)"),
    x = "UMAP 1", y = "UMAP 2", color = "TCR Type"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.6),
    axis.title = element_text(size = 13, face = "bold"),
    axis.text = element_text(size = 11),
    legend.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

# Save
ggsave(file.path(qc_folder, "Plot4bbbb_UMAP_378_Viral_TCRs_Final.png"),
       plot = umap_plot, width = 10, height = 7, dpi = 300)

# View
print(umap_plot)








# 10. Additional Plots — Clean Circular Points

# UMAP with 14 clusters using circle dots (no rasterization)
umap_plot <- DimPlot(
  combined,
  reduction = "umap",
  label = TRUE,
  raster = FALSE,  # ensures circle points
  size = 0.2       # smaller dot size for clarity
) +
  ggtitle("UMAP - 12 Samples (Exactly 14 Clusters)") +
  theme_minimal()

# Sample Contribution by Cluster (bar plot)
sample_contribution <- as.data.frame(table(Idents(combined), combined$orig.ident))
colnames(sample_contribution) <- c("Cluster", "Sample", "Count")
sample_contribution_plot <- ggplot(sample_contribution, aes(x = Cluster, y = Count, fill = Sample)) +
  geom_bar(stat = "identity", position = "fill") +
  ylab("Proportion") +
  ggtitle("Sample Contribution by Cluster (14 Clusters)") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set3")

# GZMK Expression Plot (violin)
if ("GZMK" %in% rownames(combined)) {
  gzmk_expression_plot <- VlnPlot(
    combined,
    features = "GZMK",
    group.by = "seurat_clusters",
    pt.size = 0.1
  ) +
    ggtitle("GZMK Expression in Each Cluster") +
    scale_fill_manual(values = scales::hue_pal()(14)) +
    theme_minimal()
} else {
  gzmk_expression_plot <- ggplot() +
    ggtitle("GZMK not found in dataset!") +
    theme_minimal()
  cat("\n WARNING: GZMK not found in dataset! Skipping GZMK plot...\n")
}

# Save All Plots
ggsave(file.path(qc_folder, "Plot1_Combined_UMAP_14Clusters_CIRCLES.png"),
       plot = umap_plot, width = 10, height = 6)

ggsave(file.path(qc_folder, "Plot2_Sample_Contribution_14Clusters.png"),
       plot = sample_contribution_plot, width = 10, height = 6)

ggsave(file.path(qc_folder, "Plot3_GZMK_Expression_14Clusters.png"),
       plot = gzmk_expression_plot, width = 10, height = 6)

# Completion message
cat("\n All circle-style UMAP and other plots saved successfully to", qc_folder, "\n")





# For 378

# ─────────────────────────────────────────────────────────────
# Shamim_Roshni - Unique Viral TCR UMAP Plot + Summary Report
# Date: Sys.time() will be used in output header
# ─────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(readxl)
  library(stringr)
})

# Define paths
qc_folder <- "/home/abmshamimuh/QC_Results4"
tcr_file <- "/home/abmshamimuh/Corridoni TCRs.xlsx"
seurat_file <- file.path(qc_folder, "Final_Seurat_Object_14Clusters.rds")

# Load Seurat object
combined <- readRDS(seurat_file)

# Load and clean TCR metadata
tcr_data <- read_excel(tcr_file, sheet = "counts with matches")
tcr_data <- tcr_data %>% rename_with(make.names)  # Avoid duplicate column name issues
tcr_data$barcode_clean <- str_replace(tcr_data$barcode, "-1$", "")
tcr_data$viral_status <- ifelse(is.na(tcr_data$organism), "Non-Viral", "Viral")

# Get unique viral barcode + CDR3β clone stats
unique_clones_excel <- unique(tcr_data$cdr3_beta[tcr_data$viral_status == "Viral"])
unique_barcodes_excel <- unique(tcr_data$barcode_clean[tcr_data$viral_status == "Viral"])

# Extract cleaned barcode from Seurat
combined$seurat_barcode <- substr(Cells(combined), 1, 16)  # first 16 bases = real barcode
barcode_df <- data.frame(seurat_barcode = combined$seurat_barcode, row.names = Cells(combined))

# Join with viral status
barcode_df <- left_join(barcode_df, tcr_data[, c("barcode_clean", "viral_status", "cdr3_beta")],
                        by = c("seurat_barcode" = "barcode_clean"))

barcode_df$viral_status[is.na(barcode_df$viral_status)] <- "Non-Viral"
combined$Viral <- barcode_df$viral_status
combined$CDR3beta <- barcode_df$cdr3_beta

# Subset only viral-labeled Seurat cells
viral_cells <- subset(combined, subset = Viral == "Viral")
viral_barcodes <- unique(viral_cells$seurat_barcode)

# Count unique viral barcodes and CDR3 clones in Seurat
unique_barcodes_seurat <- unique(viral_cells$seurat_barcode)
unique_clones_seurat <- unique(viral_cells$CDR3beta)
expansion <- length(Cells(viral_cells)) / length(unique_clones_seurat)

# Generate data frame for UMAP
df_all <- Embeddings(combined, "umap") %>%
  as.data.frame() %>%
  mutate(Group = ifelse(rownames(.) %in% Cells(viral_cells), "Viral TCRs", "Non-Viral TCRs")) %>%
  mutate(UMAP_1 = .[[1]], UMAP_2 = .[[2]])

# Plot UMAP
umap_plot <- ggplot(df_all, aes(x = UMAP_1, y = UMAP_2, color = Group)) +
  geom_point(size = 0.5, alpha = 0.7) +
  scale_color_manual(values = c("Viral TCRs" = "red", "Non-Viral TCRs" = "gray80")) +
  labs(title = "UMAP of Seurat Object\nViral TCR Cells in Red (n = 378)",
       x = "UMAP 1", y = "UMAP 2", color = "TCR Type") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 12)
  )

# Save UMAP
ggsave(file.path(qc_folder, "1-378Plot_Viral_TCR_UMAP_Unique_Only.png"), plot = umap_plot, width = 10, height = 7, dpi = 300)

# Generate summary
summary_txt <- c(
  paste0(" Viral TCR Summary (Generated: ", Sys.time(), ")"),
  "--------------------------------------------------",
  paste0("a) Total viral TCR entries in Excel (non-unique): ", nrow(tcr_data[tcr_data$viral_status == "Viral", ])),
  paste0("b) Total unique viral TCR clones in Excel (unique CDR3β): ", length(unique_clones_excel)),
  paste0("c) Total unique viral TCRs in Excel (by barcode): ", length(unique_barcodes_excel)),
  paste0("d) Total cells with viral unique TCRs (in Seurat object): ", length(Cells(viral_cells))),
  paste0("e) Matched unique viral TCR barcodes in Seurat object: ", length(unique_barcodes_seurat)),
  paste0("f) Unique viral TCR clones in Seurat object: ", length(unique_clones_seurat)),
  paste0("g) Average expansion per clone in Seurat object: ", round(expansion, 2))
)

# Save to file
writeLines(summary_txt, file.path(qc_folder, "1-378Viral_TCR_UMAP_Summary.txt"))

# Print preview
cat(paste(summary_txt, collapse = "\n"))




# ———————————————————————————————————————————
# Shamim_Roshni - Unique Viral TCR UMAP Plot + Summary Report
# Purpose: Highlight both total matched viral TCR cells (1134) AND
#          uniquely matched viral TCRs (378) from Excel in Seurat object
# ———————————————————————————————————————————

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(readxl)
  library(stringr)
})

# Set paths
qc_folder <- "/home/abmshamimuh/QC_Results4"
tcr_file <- "/home/abmshamimuh/Corridoni TCRs.xlsx"
seurat_file <- file.path(qc_folder, "Final_Seurat_Object_14Clusters.rds")

# Load data
combined <- readRDS(seurat_file)
tcr_data <- read_excel(tcr_file, sheet = "counts with matches") %>%
  rename_with(make.names)
tcr_data$barcode_clean <- str_replace(tcr_data$barcode, "-1$", "")
tcr_data$viral_status <- ifelse(is.na(tcr_data$organism), "Non-Viral", "Viral")

# Map Seurat barcodes
combined$seurat_barcode <- substr(Cells(combined), 1, 16)
barcode_df <- data.frame(seurat_barcode = combined$seurat_barcode, row.names = Cells(combined))
barcode_df <- left_join(barcode_df, tcr_data[, c("barcode_clean", "viral_status", "cdr3_beta")],
                        by = c("seurat_barcode" = "barcode_clean"))
barcode_df$viral_status[is.na(barcode_df$viral_status)] <- "Non-Viral"
combined$Viral <- barcode_df$viral_status
combined$CDR3beta <- barcode_df$cdr3_beta

#  Total Seurat cells with viral label (including duplicates)
all_viral_cells <- subset(combined, subset = Viral == "Viral")

#  Unique match: Keep 1 cell per unique barcode
unique_barcodes <- all_viral_cells$seurat_barcode[!duplicated(all_viral_cells$seurat_barcode)]
unique_viral_cells <- all_viral_cells[, all_viral_cells$seurat_barcode %in% unique_barcodes]

# Create UMAP DataFrame
umap_df <- Embeddings(combined, "umap") %>%
  as.data.frame() %>%
  mutate(
    UMAP_1 = .[[1]],
    UMAP_2 = .[[2]],
    Group = ifelse(rownames(.) %in% Cells(unique_viral_cells), "Viral TCRs", "Non-Viral TCRs")
  )

# Plot UMAP
umap_plot <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = Group)) +
  geom_point(size = 0.6, alpha = 0.8) +
  scale_color_manual(values = c("Viral TCRs" = "red", "Non-Viral TCRs" = "gray80")) +
  labs(
    title = "UMAP of Seurat Object\nViral TCR Cells in Red (n = 378)",
    x = "UMAP 1", y = "UMAP 2", color = "TCR Type"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 12)
  )

ggsave(file.path(qc_folder, "2-378Plot_Viral_TCR_UMAP_Unique_Only.png"), plot = umap_plot, width = 10, height = 7, dpi = 300)

# ── Summary Stats ──
total_excel_tcrs <- nrow(tcr_data[tcr_data$viral_status == "Viral", ])
unique_excel_clones <- length(unique(tcr_data$cdr3_beta[tcr_data$viral_status == "Viral"]))
unique_excel_barcodes <- length(unique(tcr_data$barcode_clean[tcr_data$viral_status == "Viral"]))
total_viral_cells <- ncol(all_viral_cells)
unique_matched_barcodes <- length(unique(unique_viral_cells$seurat_barcode))
unique_matched_clones <- length(unique(unique_viral_cells$CDR3beta))
avg_expansion <- round(total_viral_cells / unique_matched_clones, 2)

summary_txt <- c(
  paste0("\U0001F4CA Viral TCR Summary (Generated: ", Sys.time(), ")"),
  "--------------------------------------------------",
  paste0("a) Total viral TCR entries in Excel (non-unique): ", total_excel_tcrs),
  paste0("b) Total unique viral TCR clones in Excel (unique CDR3β): ", unique_excel_clones),
  paste0("c) Total unique viral TCRs in Excel (by barcode): ", unique_excel_barcodes),
  paste0("d) Total Seurat cells matched to viral barcodes in Excel (all viral-labeled): ", total_viral_cells),
  paste0("e) Matched unique viral TCR barcodes in Seurat object: ", unique_matched_barcodes),
  paste0("f) Unique viral TCR clones in Seurat object: ", unique_matched_clones),
  paste0("g) Average expansion per clone in Seurat object: ", avg_expansion)
)

writeLines(summary_txt, file.path(qc_folder, "2-378Viral_TCR_UMAP_Summary.txt"))
cat(paste(summary_txt, collapse = "\n"))










# Below are okay/no error


# Required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(stringr)
})

# Set output directory
output_dir <- file.path("/home/abmshamimuh/QC_Results4")

# ️⃣ Cluster UMAP with Viral Overlay (fully fixed)

# Base cluster UMAP
umap_cluster <- DimPlot(combined, reduction = "umap", group.by = "seurat_clusters", label = TRUE,
                        raster = FALSE, size = 0.2) +
  ggtitle("UMAP (Clusters) with Viral TCR Overlay") +
  theme_minimal()

#  Define viral cells
viral_cells <- subset(combined, subset = Viral == "Viral")

#  Extract UMAP and rename columns
df_viral <- Embeddings(viral_cells, reduction = "umap") %>%
  as.data.frame() %>%
  rename(UMAP_1 = umap_1, UMAP_2 = umap_2)

#  Overlay viral red dots
umap_cluster_overlay <- umap_cluster +
  geom_point(data = df_viral, aes(x = UMAP_1, y = UMAP_2), color = "red", size = 0.3, inherit.aes = FALSE)

# Save final image
ggsave(file.path(output_dir, "Plot1A_UMAP_Clusters_With_ViralOverlay.png"),
       plot = umap_cluster_overlay, width = 10, height = 7)

#Viral TCRs by Cluster (Bar Plot)
viral_cluster_dist <- table(combined$seurat_clusters[combined$Viral == "Viral"])
viral_cluster_df <- as.data.frame(viral_cluster_dist)
colnames(viral_cluster_df) <- c("Cluster", "Viral_Count")

viral_bar <- ggplot(viral_cluster_df, aes(x = Cluster, y = Viral_Count, fill = Cluster)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  ggtitle("Viral TCR Count per Cluster") +
  guides(fill = "none")

ggsave(file.path(output_dir, "Bar_ViralTCRs_by_Cluster.png"), viral_bar, width = 8, height = 6)

#  Viral TCR Pie/Bar by Sample (and HC/UC if available)
viral_sample_dist <- table(combined$orig.ident[combined$Viral == "Viral"])
viral_sample_df <- as.data.frame(viral_sample_dist)
colnames(viral_sample_df) <- c("Sample", "Viral_Count")

viral_pie <- ggplot(viral_sample_df, aes(x = "", y = Viral_Count, fill = Sample)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  ggtitle("Viral TCR Proportion by Sample") +
  theme_void()

viral_bar2 <- ggplot(viral_sample_df, aes(x = Sample, y = Viral_Count, fill = Sample)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  ggtitle("Viral TCR Count by Sample")

ggsave(file.path(output_dir, "Pie_ViralTCRs_by_Sample.png"), viral_pie, width = 6, height = 6)
ggsave(file.path(output_dir, "Bar_ViralTCRs_by_Sample.png"), viral_bar2, width = 8, height = 6)


#For percentage

# Calculate counts and percentages
viral_cluster_df <- table(combined$seurat_clusters[combined$Viral == "Viral"]) %>%
  as.data.frame()
colnames(viral_cluster_df) <- c("Cluster", "Viral_Count")
viral_cluster_df$Percentage <- 100 * viral_cluster_df$Viral_Count / total_viral_cells

viral_bar <- ggplot(viral_cluster_df, aes(x = Cluster, y = Percentage, fill = Cluster)) +
  geom_bar(stat = "identity", color = "black") +
  theme_minimal(base_size = 14) +
  labs(
    title = paste0("Viral TCRs per Cluster (Percentage, n = ", total_viral_cells, ")"),
    x = "Cluster", y = "Percentage of Viral TCRs"
  ) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "none"
  )

ggsave(file.path(output_dir, "Plot5_Bar_ViralTCRs_by_Cluster_Percent.png"), viral_bar, width = 8, height = 6)

# Calculate counts and percentages
viral_sample_df <- table(combined$orig.ident[combined$Viral == "Viral"]) %>%
  as.data.frame()
colnames(viral_sample_df) <- c("Sample", "Viral_Count")
viral_sample_df$Percentage <- 100 * viral_sample_df$Viral_Count / total_viral_cells
viral_sample_df$Label <- paste0(round(viral_sample_df$Percentage, 1), "%")

viral_pie <- ggplot(viral_sample_df, aes(x = "", y = Percentage, fill = Sample)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y") +
  geom_text(aes(label = Label), position = position_stack(vjust = 0.5), size = 4, fontface = "bold") +
  labs(title = paste0("Viral TCR Distribution by Sample (Pie, n = ", total_viral_cells, ")")) +
  theme_void(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.title = element_blank()
  )

ggsave(file.path(output_dir, "Plot7_Pie_ViralTCRs_by_Sample_Percent.png"), viral_pie, width = 7, height = 6)

# Bar Plot
viral_bar2 <- ggplot(viral_sample_df, aes(x = Sample, y = Percentage, fill = Sample)) +
  geom_bar(stat = "identity", color = "black") +
  labs(
    title = paste0("Viral TCRs per Sample (Percentage, n = ", total_viral_cells, ")"),
    x = "Sample", y = "Percentage of Viral TCRs"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "none"
  )

ggsave(file.path(output_dir, "Plot6_Bar_ViralTCRs_by_Sample_Percent.png"), viral_bar2, width = 9, height = 6)






# ────────────────────────────────────────────────
#  Plot 9: UMAP Overlay of GZMK Expression
# ────────────────────────────────────────────────

gzmk_umap <- FeaturePlot(
  combined,
  features = "GZMK",
  cols = c("lightgray", "darkgreen"),
  pt.size = 0.5
) +
  ggtitle("UMAP Overlay: GZMK Expression") +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 11)
  )

ggsave(
  filename = file.path(output_dir, "Plot9_GZMK_UMAP.png"),
  plot = gzmk_umap, width = 8, height = 6, dpi = 300
)

# ────────────────────────────────────────────────
#  Plot 10b: Violin Plot of GZMK by Viral Status
# ────────────────────────────────────────────────

# Step 1: Fetch GZMK expression & viral status
gzmk_data <- FetchData(combined, vars = c("GZMK", "Viral")) %>% na.omit()

# Step 2: Wilcoxon rank-sum test
gzmk_wilcox <- wilcox.test(GZMK ~ Viral, data = gzmk_data)

# Step 3: P-value label
p_val <- gzmk_wilcox$p.value
p_label <- if (p_val < 0.001) "p < 0.001" else paste0("p = ", signif(p_val, 3))

# Step 4: Violin plot with annotation
violin_plot <- VlnPlot(
  combined,
  features = "GZMK",
  group.by = "Viral",
  pt.size = 0.2,
  cols = c("gray70", "red3")
) +
  ggtitle("GZMK Expression in Viral vs Non-Viral TCR Cells") +
  ylab("GZMK Expression") + xlab("TCR Type") +
  annotate("text",
           x = 1.5,
           y = max(gzmk_data$GZMK, na.rm = TRUE) * 1.05,
           label = p_label,
           size = 5,
           fontface = "italic") +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold")
  )

# Step 5: Save violin plot
ggsave(
  filename = file.path(output_dir, "Plot10b_GZMK_Viral_Violin_Annotated.png"),
  plot = violin_plot, width = 7, height = 6, dpi = 300
)
# this is completely expected in scaled or SCTransformed data. 
# It doesn't mean expression is biologically negative — it's just a statistical normalization artifact.





# to visualize raw log-normalized GZMK expression (all ≥ 0), you can:
# Use RNA assay to avoid scaled/negative values
DefaultAssay(combined) <- "RNA"

# Step 1: Extract raw expression and metadata
gzmk_data <- FetchData(combined, vars = c("GZMK", "Viral")) %>%
  na.omit()

# Step 2: Run Wilcoxon rank-sum test
gzmk_wilcox <- wilcox.test(GZMK ~ Viral, data = gzmk_data)

# Step 3: Format p-value label
p_val <- gzmk_wilcox$p.value
p_label <- if (p_val < 0.001) {
  "p < 0.001"
} else {
  paste0("p = ", signif(p_val, 3))
}

# Step 4: Violin Plot
violin_plot <- VlnPlot(
  combined,
  features = "GZMK",
  group.by = "Viral",
  pt.size = 0.1,
  assay = "RNA"
) +
  ggtitle("GZMK Expression in Viral vs Non-Viral TCR Cells") +
  xlab("TCR Type") +
  ylab("GZMK Expression") +
  theme_classic(base_size = 16) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  ) +
  annotate(
    "text",
    x = 1.5,
    y = max(gzmk_data$GZMK, na.rm = TRUE) * 1.05,
    label = p_label,
    size = 5,
    fontface = "italic"
  )

# Step 5: Save
ggsave(
  file.path(output_dir, "Plot10bb_GZMK_Viral_Violin_Annotated.png"),
  plot = violin_plot,
  width = 8, height = 6, dpi = 300
)


library(dplyr)
library(ggplot2)

# 1. Summarize counts for percentage labels
viral_counts <- gzmk_data %>%
  group_by(Viral) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(Percentage = round(n / sum(n) * 100, 1),
         Label = paste0(Viral, " (", Percentage, "%)"))

# 2. Match labels back to main data
gzmk_data <- gzmk_data %>%
  left_join(viral_counts[, c("Viral", "Label")], by = "Viral")

# 3. Plot
p <- ggplot(gzmk_data, aes(x = Label, y = GZMK, fill = Label)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6) +
  geom_jitter(width = 0.2, size = 0.5, alpha = 0.3) +
  scale_fill_manual(values = c(
    "Non-Viral (x%)" = "gray",  # These will be overwritten by actual labels below
    "Viral (x%)" = "red"
  )) +
  labs(
    title = "GZMK Expression by TCR Type",
    x = "TCR Type", y = "GZMK Expression", fill = "TCR Status"
  ) +
  theme_classic(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    legend.position = "right",
    legend.title = element_text(face = "bold")
  )

# 4. Save
ggsave(file.path(output_dir, "Plot11a_Boxplot_GZMK_Viral_PercentLegend.png"),
       plot = p, width = 8, height = 6, dpi = 300)

print(p)


library(dplyr)
library(ggplot2)
library(scales)

# Extract GZMK expression vector
#gzmk_expr <- GetAssayData(combined, slot = "data")["GZMK", ]
# Detect assay
assay_name <- if ("SCT" %in% Assays(combined)) "SCT" else DefaultAssay(combined)

# Detect available layers
available_layers <- Layers(combined[[assay_name]])

# Choose best available layer
target_layer <- if ("data" %in% available_layers) {
  "data"
} else if ("counts" %in% available_layers) {
  "counts"
} else {
  stop(" No usable layer found in assay: ", assay_name)
}

# Extract GZMK expression
gzmk_expr <- GetAssayData(combined, assay = assay_name, layer = target_layer)["GZMK", ]
combined$GZMK_expr <- as.numeric(gzmk_expr[Cells(combined)])

# Prepare data for plotting
gzmk_data <- FetchData(combined, vars = c("GZMK_expr", "Viral")) %>%
  rename(GZMK = GZMK_expr) %>%
  na.omit()

# Add to metadata
combined$GZMK_Pos <- ifelse(gzmk_expr > 0, 1, 0)

# Compute cluster-wise percentages
gzmk_cluster_df <- combined@meta.data %>%
  group_by(Cluster = seurat_clusters) %>%
  summarise(
    Total = n(),
    GZMK_Positive = sum(GZMK_Pos),
    Percent = round(100 * GZMK_Positive / Total, 1),
    .groups = "drop"
  )

# Plot
plot_12a <- ggplot(gzmk_cluster_df, aes(x = Cluster, y = Percent, fill = Cluster)) +
  geom_bar(stat = "identity", color = "black") +
  labs(
    title = "GZMK+ Cell Contribution by Cluster",
    x = "Cluster", y = "Percentage of GZMK+ Cells"
  ) +
  scale_y_continuous(labels = percent_format(scale = 1)) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "none"
  )

ggsave(file.path(output_dir, "Plot12a_GZMK_Contribution_by_Cluster.png"),
       plot = plot_12a, width = 8, height = 6, dpi = 300)



# statistical test
# Load required packages
library(dplyr)
library(Seurat)
library(ggplot2)
library(scales)

# Step 1: Extract GZMK expression and add metadata
assay_name <- if ("SCT" %in% Assays(combined)) "SCT" else DefaultAssay(combined)
available_layers <- Layers(combined[[assay_name]])

target_layer <- if ("data" %in% available_layers) {
  "data"
} else if ("counts" %in% available_layers) {
  "counts"
} else {
  stop("No usable layer found.")
}

gzmk_expr <- GetAssayData(combined, assay = assay_name, layer = target_layer)["GZMK", ]
combined$GZMK_expr <- as.numeric(gzmk_expr[Cells(combined)])
combined$GZMK_Pos <- ifelse(combined$GZMK_expr > 0, 1, 0)

# Step 2: Summary of GZMK+ percentage by cluster
gzmk_cluster_summary <- combined@meta.data %>%
  group_by(Cluster = seurat_clusters) %>%
  summarise(
    Total_Cells = n(),
    GZMK_Positive = sum(GZMK_Pos),
    Percent_Positive = round(100 * GZMK_Positive / Total_Cells, 1),
    .groups = "drop"
  )
print(gzmk_cluster_summary)

# Step 3a: Fisher’s Exact Test between Cluster 4 and Cluster 6
tbl_fisher <- table(combined$seurat_clusters, combined$GZMK_Pos)
fisher_pair <- fisher.test(tbl_fisher[c("4", "6"), ])
cat("\n Plot12a_Fisher’s Exact Test (Cluster 4 vs 6): p-value =", signif(fisher_pair$p.value, 4), "\n")

# Step 3b: Chi-Square Test across all 14 clusters
chisq_result <- chisq.test(tbl_fisher)
cat("Plot12a_Chi-Square Test across all clusters (0–13): p-value =", signif(chisq_result$p.value, 4), "\n")

# Format and print p-value with higher precision
formatted_pval <- format.pval(chisq_result$p.value, digits = 10, scientific = TRUE)
cat("Plot12a_Formatted p-value (scientific notation):", formatted_pval, "\n")

# For save the results

# Step 3a: Fisher’s Exact Test between Cluster 4 and Cluster 6
tbl_fisher <- table(combined$seurat_clusters, combined$GZMK_Pos)
fisher_pair <- fisher.test(tbl_fisher[c("4", "6"), ])

# Step 3b: Chi-Square Test across all 14 clusters
chisq_result <- chisq.test(tbl_fisher)

# Format Chi-Square p-value with higher precision
formatted_pval <- format.pval(chisq_result$p.value, digits = 10, scientific = TRUE)

# Step 4: Save results to a text file
result_text <- paste0(
  "Statistical Test Results for GZMK+ by Cluster (Plot12a)\n",
  "--------------------------------------------------------\n",
  "Fisher’s Exact Test (Cluster 4 vs 6):\n",
  "  p-value = ", signif(fisher_pair$p.value, 4), "\n\n",
  "Chi-Square Test across all clusters (0–13):\n",
  "  Raw p-value = ", signif(chisq_result$p.value, 4), "\n",
  "  Scientific notation = ", formatted_pval, "\n"
)

# Save to file
writeLines(result_text, "/home/abmshamimuh/QC_Results4/Plot12a_StatisticalTests.txt")

# Also print to screen
cat(result_text)










# Count GZMK+ and total cells per sample
gzmk_sample_df <- combined@meta.data %>%
  group_by(Sample = orig.ident) %>%
  summarise(
    Total = n(),
    GZMK_Positive = sum(GZMK_Pos),
    Percent = round(100 * GZMK_Positive / Total, 1),
    .groups = "drop"
  )

plot_12b <- ggplot(gzmk_sample_df, aes(x = Sample, y = Percent, fill = Sample)) +
  geom_bar(stat = "identity", color = "black") +
  labs(
    title = "GZMK+ Cell Contribution by Sample",
    x = "Sample", y = "Percentage of GZMK+ Cells"
  ) +
  scale_y_continuous(labels = percent_format(scale = 1)) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

ggsave(file.path(output_dir, "Plot12b_GZMK_Contribution_by_Sample.png"),
       plot = plot_12b, width = 10, height = 6, dpi = 300)




# UMAP Colored by Sample
umap_sample <- DimPlot(combined, reduction = "umap", group.by = "orig.ident",
                       raster = FALSE, size = 0.2) +
  ggtitle("UMAP Colored by Sample") +
  theme_minimal()

ggsave(file.path(output_dir, "Plot1B_UMAP_Colored_by_Sample.png"), umap_sample, width = 10, height = 7)










#  Setup and Library Loading
options(future.globals.maxSize = 450 * 1024^3)
options(repos = c(CRAN = "https://cran.r-project.org"))
.libPaths(c("/home/abmshamimuh/R/x86_64-pc-linux-gnu-library/4.4", .libPaths()))

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(readxl)
  library(stringr)
  library(RColorBrewer)
  library(scales)
  library(tidyr)
  library(tibble)
  library(lme4)
  library(emmeans)
  library(ggpubr)
  library(eulerr)
})

# Data Loading

# Paths
data_folder <- "/home/abmshamimuh/"
qc_folder <- file.path(data_folder, "QC_Results4")
tcr_file <- file.path(data_folder, "Corridoni TCRs.xlsx")
seurat_rds <- file.path(qc_folder, "Final_Seurat_Object_14Clusters.rds")
output_dir <- qc_folder

# Load Seurat object
combined <- readRDS(seurat_rds)
combined$Cluster14 <- combined$seurat_clusters  # Ensure Cluster14 is available

# Load and clean TCR metadata
tcr_data <- read_excel(tcr_file, sheet = "counts with matches") %>% 
  rename_with(make.names)
tcr_data$barcode_clean <- str_replace(tcr_data$barcode, "-1$", "")
tcr_data$viral_status <- ifelse(is.na(tcr_data$organism), "Non-Viral", "Viral")

# Extract barcode and annotate combined Seurat object
combined$seurat_barcode <- substr(Cells(combined), 1, 16)
barcode_df <- data.frame(seurat_barcode = combined$seurat_barcode, row.names = Cells(combined))
barcode_df <- left_join(barcode_df, tcr_data[, c("barcode_clean", "viral_status", "cdr3_beta", "Sample", "organism")],
                        by = c("seurat_barcode" = "barcode_clean"))
barcode_df$viral_status[is.na(barcode_df$viral_status)] <- "Non-Viral"
combined$Viral <- barcode_df$viral_status
combined$cdr3_beta <- barcode_df$cdr3_beta
combined$organism <- barcode_df$organism
combined$sample <- combined$orig.ident  # Add sample column for compatibility

#  Plot 22: GZMK Expression by Cluster
plot_22 <- VlnPlot(combined, features = "GZMK", group.by = "Cluster14", pt.size = 0.1) +
  ggtitle("GZMK Expression by Cluster") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave(file.path(output_dir, "Plot22_GZMK_Expression_by_Cluster.png"), plot = plot_22, width = 10, height = 6)

#  Plot 33: Sample Contribution by Cluster
plot_33 <- ggplot(combined@meta.data, aes(x = Cluster14, fill = orig.ident)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = percent) +
  labs(title = "Sample Contribution by Cluster", y = "Percentage", x = "Cluster") +
  theme_minimal(base_size = 14) +
  theme(axis.title = element_text(face = "bold"), legend.title = element_blank())
ggsave(file.path(output_dir, "Plot33_Sample_Contribution_by_Cluster.png"), plot = plot_33, width = 9, height = 6)

#  Plot 44: Viral vs Non-Viral Composition by Cluster
plot_44 <- ggplot(combined@meta.data, aes(x = Cluster14, fill = Viral)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = c("Non-Viral" = "gray", "Viral" = "red")) +
  scale_y_continuous(labels = percent) +
  labs(title = "Viral vs Non-Viral Composition", y = "Percentage", x = "Cluster") +
  theme_minimal(base_size = 14) +
  theme(axis.title = element_text(face = "bold"), legend.title = element_blank())
ggsave(file.path(output_dir, "Plot44_Viral_vs_NonViral_Cluster_Composition.png"), plot = plot_44, width = 9, height = 6)

#  Plot 44 Statistics (Chi-square test)
plot4_stats <- combined@meta.data %>%
  group_by(Cluster14) %>%
  summarise(
    viral_cells = sum(Viral == "Viral"),
    total_cells = n(),
    p_value = chisq.test(table(Viral))$p.value
  ) %>%
  mutate(
    fdr = p.adjust(p_value, method = "BH"),
    significance = case_when(
      fdr < 0.001 ~ "***",
      fdr < 0.01 ~ "**",
      fdr < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )
write.csv(plot4_stats, file.path(output_dir, "Stats_Plot44_Cluster_Viral_Composition.csv"), row.names = FALSE)





