# ---- Setup for RCC ----
options(future.globals.maxSize = 450 * 1024^3)

required_packages <- c("Seurat", "ggplot2", "dplyr", "patchwork", "readxl", "stringr")
invisible(lapply(required_packages, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg, repos = "https://cran.r-project.org")
}))
suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(dplyr); library(patchwork)
  library(readxl); library(stringr)
})

# ---- Paths ----
input_dir <- "/home/abmshamimuh/QC_Results4"
tcr_file <- file.path(input_dir, "Filtered_TCR_Data__S21_S34_.csv")

# ---- Load TCR Data ----
tcr_data <- read.csv(tcr_file)

# ---- Sample IDs and Seurat Files ----
sample_ids <- c(
  "GSM4339769_C141", "GSM4339770_C142", "GSM4339771_C143",
  "GSM4339772_C144", "GSM4339773_C145", "GSM4339774_C146",
  "GSM4475048_C51",  "GSM4475049_C52",  "GSM4475050_C100",
  "GSM4475051_C148", "GSM4475052_C149", "GSM4475053_C152"
)
seurat_files <- file.path(input_dir, paste0(sample_ids, "_seurat_final.rds"))

# ---- Load and Merge Seurat Objects ----
cat("📦 Loading Seurat samples...\n")
seurat_list <- lapply(seurat_files, readRDS)
cat("🔁 Merging Seurat objects...\n")
combined <- merge(seurat_list[[1]], y = seurat_list[-1])

# ---- Step 1: Unique Viral TCRs in Excel ----
unique_viral_tcrs_excel <- unique(tcr_data$barcode)
cat("1️⃣ Unique Viral TCRs in Excel:", length(unique_viral_tcrs_excel), "\n")

# ---- Step 2: Match TCR Barcodes to Seurat ----
tcr_data$barcode_clean <- tcr_data$barcode
tcr_data$combined_key <- NA_character_

for (sample_id in unique(tcr_data$Sample)) {
  matched_sample <- grep(sample_id, sample_ids, value = TRUE)
  if (length(matched_sample) == 1) {
    sel <- tcr_data$Sample == sample_id
    tcr_data$combined_key[sel] <- paste0(tcr_data$barcode_clean[sel], "-1_", matched_sample)
  }
}

matched_keys <- intersect(tcr_data$combined_key, colnames(combined))
unique_viral_matched <- length(unique(matched_keys))
cat("2️⃣ Unique Viral TCRs matched to Seurat:", unique_viral_matched, "\n")

# ---- Step 3: Clonal Viral TCRs in Seurat ----
tcr_matches <- tcr_data %>%
  filter(combined_key %in% matched_keys) %>%
  select(combined_key, raw_clonotype_id)

clonal_counts <- tcr_matches %>%
  group_by(raw_clonotype_id) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  filter(n_cells > 1)

cat("3️⃣ Clonal Viral TCRs in Seurat (>1 cell):", nrow(clonal_counts), "\n")

# ---- Step 4: Unique Clonal Seurat Barcodes ----
clonal_ids <- clonal_counts$raw_clonotype_id

unique_clonal_barcodes <- tcr_matches %>%
  filter(raw_clonotype_id %in% clonal_ids) %>%
  pull(combined_key) %>%
  unique()

cat("4️⃣ Unique Viral TCRs (clonal) matched to Seurat:", length(unique_clonal_barcodes), "\n")

# ---- Step 5: Total Unique Viral Barcodes Matched (clonal + single) ----
total_matched_unique_barcodes <- unique(tcr_matches$combined_key)
cat("5️⃣ Unique Viral TCRs (clonal + single) matched to Seurat:", length(total_matched_unique_barcodes), "\n")

# ---- Step 6: Total Matched Seurat Cells (duplicates included) ----
# Note: If any barcode appears in multiple Seurat cells, it will be counted for each.
viral_seurat_cells <- tcr_data %>%
  filter(combined_key %in% colnames(combined)) %>%
  pull(combined_key)

cat("6️⃣ Total Viral TCR-Labeled Seurat Cells (including duplicates):", length(viral_seurat_cells), "\n")
