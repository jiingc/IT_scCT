# Signac-0.R: Load and process each replicate independently for each histone mark
# This script constructs Seurat objects for each replicate, then merges them within each mark.

suppressPackageStartupMessages({
  library(Signac)
  library(Seurat)
  library(GenomeInfoDb)
  library(EnsDb.Mmusculus.v79)
  library(GenomicRanges)
  library(Matrix)
})

### Set parameters
PROJECT <- "ITscCT"
TISSUE <- "MG"
TIMEPOINT <- "10W"
GENOME <- "mm10"
BASE_DIR <- "/path/to/project"
PEAK_DIR <- file.path(BASE_DIR, "peaks")
FRAG_DIR <- file.path(BASE_DIR, "fragments")
REPLICATES <- c("Rep1", "Rep2")

### Helper to load one replicate for one mark
load_scCutTag_replicate <- function(mark, replicate) {
  sample_id <- paste(PROJECT, TISSUE, TIMEPOINT, replicate, mark, sep = "-")
  frag_file <- file.path(FRAG_DIR, replicate, paste0(sample_id, ".fragments.sort.bed.gz"))
  peak_file <- file.path(PEAK_DIR, paste0(sample_id, "_peaks.broadPeak"))

  # Load peak and fragments
  peaks_df <- read.table(
    file = peak_file,
    col.names = c("chr", "start", "end", "name", "score", "strand", "fc", "-log10p", "-log10q")
  )
  peaks_gr <- makeGRangesFromDataFrame(peaks_df, keep.extra.columns = TRUE)
  frag_obj <- CreateFragmentObject(frag_file, genome = GENOME)

  # Count matrix and assay
  mat <- FeatureMatrix(fragments = frag_obj, features = peaks_gr)
  mat <- mat[rowSums(mat) > 0, ]
  assay <- CreateChromatinAssay(counts = mat, fragments = frag_obj, genome = GENOME)

  # Seurat object
  seu <- CreateSeuratObject(counts = assay, assay = "peaks")
  seu$replicate <- replicate
  seu$mark <- mark
  seu$sample_id <- sample_id

  # Add FRiP and fragment count
  total_frags <- CountFragments(frag_file)
  rownames(total_frags) <- total_frags$CB
  seu$fragments <- total_frags[colnames(seu), "frequency_count"]
  seu <- FRiP(seu, assay = "peaks", total.fragments = "fragments", col.name = "FRiP")

  return(seu)
}

### Load all replicates for each mark
marks <- c("H3K4me3", "H3K27ac")
seurat_list <- list()

for (mark in marks) {
  rep_list <- lapply(REPLICATES, function(rep) load_scCutTag_replicate(mark, rep))
  merged <- merge(rep_list[[1]], y = rep_list[-1])
  seurat_list[[mark]] <- merged
}

### Gene annotations
annotations <- GetGRangesFromEnsDb(EnsDb.Mmusculus.v79)
seqlevels(annotations) <- paste0("chr", seqlevels(annotations))
genome(annotations) <- GENOME

for (mark in marks) {
  Annotation(seurat_list[[mark]]) <- annotations
}

### Filtering step
for (mark in marks) {
  obj <- seurat_list[[mark]]
  filtered_cells <- rownames(obj@meta.data[obj$fragments > 50 & obj$FRiP > 0.3, ])
  seurat_list[[mark]] <- subset(obj, cells = filtered_cells)
}

### Output objects
H3K4me3_seurat_filtered <- seurat_list[["H3K4me3"]]
H3K27ac_seurat_filtered <- seurat_list[["H3K27ac"]]

saveRDS(H3K4me3_seurat_filtered, file = file.path(BASE_DIR, "rds/H3K4me3_seurat_filtered.rds"))
saveRDS(H3K27ac_seurat_filtered, file = file.path(BASE_DIR, "rds/H3K27ac_seurat_filtered.rds"))

message("All histone mark Seurat objects loaded, merged and filtered.")

