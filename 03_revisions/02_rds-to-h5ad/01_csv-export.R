# export_seurat_layers_to_csv.R

library(Seurat)
library(SeuratObject)

# Set your source and output directories
input_dir <- "/path/to/soupx_rds"
output_base <- "/path/to/soupx_export"

samples <- c("1", "3", "5", "6", "7", "8")

dir.create(output_base, showWarnings = FALSE)

for (sample in samples) {
  message("Processing sample ", sample)

  input_file <- file.path(input_dir, paste0(sample, "_soupx_corrected.rds"))
  output_dir <- file.path(output_base, sample)
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  sobj <- readRDS(input_file)

  # Extract and write SoupX-corrected counts (transpose for Python compatibility)
  corrected <- GetAssayData(sobj, layer = "counts")
  write.csv(t(as.matrix(corrected)), file = file.path(output_dir, "counts_corrected_soupx.csv"))

  # Extract and write original (pre-SoupX) counts (transpose for Python compatibility)
  original <- GetAssayData(sobj, layer = "original")
  write.csv(t(as.matrix(original)), file = file.path(output_dir, "counts_raw.csv"))

  # Metadata
  write.csv(sobj@meta.data, file = file.path(output_dir, "metadata.csv"))

  # Features (genes)
  features <- data.frame(gene = rownames(corrected))
  write.csv(features, file = file.path(output_dir, "features.csv"), row.names = FALSE)

  message("Exported: ", output_dir)
}

# Log any warnings encountered
if (length(warnings()) > 0) {
  message("\n========== WARNINGS ==========")
  print(warnings())
  message("================================")
}