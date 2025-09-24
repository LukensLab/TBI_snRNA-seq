library(Seurat)
library(SoupX)
library(SeuratObject)

setwd("/path/to/cellranger")

mad_outlier <- function(sobj, metric, nmads){
  M <- sobj@meta.data[[metric]]
  median_M <- median(M, na.rm = TRUE)
  mad_M <- mad(M, na.rm = TRUE)
  outlier <- (M < (median_M - nmads * mad_M)) | (M > (median_M + nmads * mad_M))
  return(outlier)
}

pp <- function(sample_id){
  path <- paste0(sample_id, "/filtered_feature_bc_matrix.h5")
  sobj <- Read10X_h5(filename = path)
  sobj <- CreateSeuratObject(counts = sobj, min.cells = 0, min.features = 200)
  sobj$sample_id <- sample_id

  sobj$log1p_total_counts <- log1p(sobj@meta.data$nCount_RNA)
  sobj$log1p_n_genes_by_counts <- log1p(sobj@meta.data$nFeature_RNA)
  sobj[["percent.mt"]] <- PercentageFeatureSet(sobj, pattern = "^mt-")

  bool_vector <- !mad_outlier(sobj, 'log1p_total_counts', 5) &
                 !mad_outlier(sobj, 'log1p_n_genes_by_counts', 5) &
                 !mad_outlier(sobj, 'percent.mt', 3)
  sobj <- subset(sobj, cells = which(bool_vector))
  return(sobj)
}

samples <- c("1", "3", "5", "6", "7", "8")
data_list <- lapply(samples, pp)
names(data_list) <- samples

get_soup_groups <- function(sobj){
  sobj <- NormalizeData(sobj, verbose = FALSE)
  sobj <- FindVariableFeatures(sobj, nfeatures = 2000, selection.method = 'vst', verbose = FALSE)
  sobj <- ScaleData(sobj, verbose = FALSE)
  sobj <- RunPCA(sobj, npcs = 20, verbose = FALSE)
  sobj <- FindNeighbors(sobj, dims = 1:20, verbose = FALSE)
  sobj <- FindClusters(sobj, resolution = 0.5, verbose = FALSE)
  return(sobj@meta.data[['seurat_clusters']])
}

add_soup_groups <- function(sobj){
  sobj$soup_group <- get_soup_groups(sobj)
  return(sobj)
}

data_list <- lapply(data_list, add_soup_groups)

make_soup <- function(sobj){
  sample_id <- as.character(sobj$sample_id[1])
  path <- paste0(sample_id, "/raw_feature_bc_matrix.h5")
  raw <- Read10X_h5(filename = path)

  toc <- GetAssayData(sobj, layer = "counts")
  sc <- SoupChannel(tod = raw, toc = toc)
  sc <- setClusters(sc, sobj$soup_group)
  sc <- autoEstCont(sc, doPlot = FALSE)
  out <- adjustCounts(sc, roundToInt = TRUE)

  # Safely set layers using Seurat v5 API
  DefaultAssay(sobj) <- "RNA"
  sobj <- SetAssayData(sobj, layer = "original", new.data = toc)
  sobj <- SetAssayData(sobj, layer = "counts", new.data = out)

  return(sobj)
}

data_list <- lapply(data_list, make_soup)

output_dir <- "/path/to/soupx_rds"
dir.create(output_dir, showWarnings = FALSE)

for (sample in names(data_list)) {
  sobj <- data_list[[sample]]

  # Correction ratio
  original <- sum(GetAssayData(sobj, layer = "original"))
  corrected <- sum(GetAssayData(sobj, layer = "counts"))
  ratio <- round(corrected / original, 3)
  message("Correction ratio (sample ", sample, "): ", ratio)

  # Save
  save_path <- file.path(output_dir, paste0(sample, "_soupx_corrected.rds"))
  saveRDS(sobj, file = save_path)
  message("Saved: ", save_path)
}