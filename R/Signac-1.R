# Signac-1.R: Add gene activity matrix, marker scoring, and dimensional reduction

suppressPackageStartupMessages({
  library(Signac)
  library(Seurat)
  library(GenomeInfoDb)
  library(EnsDb.Mmusculus.v79)
  library(GenomicRanges)
  library(Matrix)
  library(dplyr)
})

#### Parameters ####
PROJECT <- "ITscCT"
TISSUE <- "MG"
TIMEPOINT <- "10W"
GENOME <- "mm10"
BASE_DIR <- "/path/to/project"
MARKS <- c("H3K4me3", "H3K27ac")
MARKER_FILE <- "../merge_Rproject/cell_marker_feature.txt"

#### Load gene annotations and extend gene bodies ####
gene.coords <- ensembldb::genes(EnsDb.Mmusculus.v79, filter = ~ gene_biotype == "protein_coding")
seqlevelsStyle(gene.coords) <- 'UCSC'
genebody.coords <- keepStandardChromosomes(gene.coords, pruning.mode = 'coarse')
genebody.coords.flat <- GenomicRanges::reduce(x = genebody.coords)
genebodyandpromoter.coords.flat <- Signac::Extend(genebody.coords.flat, upstream = 10000)
genebodyandpromoter.coords.flat$name <- gene.coords[nearest(genebodyandpromoter.coords.flat, genebody.coords)]$gene_name

#### Load marker gene list ####
features_all <- readLines(MARKER_FILE)

#### Define marker sets ####
features_Mature_luminal <- c("Krt18","Krt8","Prlr", "Foxa1", "Cited1")
features_Luminal_progenitor <- c("Aldh1a3","Cd14","Elf5","Kit")
features_combined <- unique(c(features_Mature_luminal, features_Luminal_progenitor))

#### Process each mark ####
for (mark in MARKS) {
  seurat_obj_name <- paste0(mark, "_seurat_filtered")
  seurat_obj <- get(seurat_obj_name)

  # Use peak assay first for scaling
  seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj), assay = "peaks")

  # Load corresponding fragment file
  sample_id <- paste(PROJECT, TISSUE, TIMEPOINT, "Rep1", mark, sep = "-")
  frag_file <- file.path(BASE_DIR, "fragments", "Rep1", paste0(sample_id, ".fragments.sort.bed.gz"))
  frag_obj <- CreateFragmentObject(frag_file, genome = GENOME)

  # Gene activity matrix
  gene.activities <- FeatureMatrix(
    fragments = frag_obj,
    features = genebodyandpromoter.coords.flat
  )
  cells_use <- intersect(colnames(seurat_obj), colnames(gene.activities))
  gene.activities <- gene.activities[, cells_use]

  gene.key <- genebodyandpromoter.coords.flat$name
  names(gene.key) <- GRangesToString(grange = genebodyandpromoter.coords.flat)
  rownames(gene.activities) <- make.unique(gene.key[rownames(gene.activities)])
  gene.activities <- gene.activities[rownames(gene.activities) != "", ]
  gene.activities <- gene.activities[rowSums(gene.activities) > 0, ]

  seurat_obj[['RNA']] <- CreateAssayObject(counts = gene.activities)
  seurat_obj <- NormalizeData(seurat_obj, assay = 'RNA', normalization.method = 'LogNormalize',
                              scale.factor = median(seurat_obj$nCount_RNA))
  DefaultAssay(seurat_obj) <- 'RNA'
  seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj), assay = "RNA")

  # Add marker scores
  f_all <- intersect(features_all, rownames(seurat_obj))
  f_lum <- intersect(features_Mature_luminal, rownames(seurat_obj))
  f_prog <- intersect(features_Luminal_progenitor, rownames(seurat_obj))
  f_combined <- intersect(features_combined, rownames(seurat_obj))

  seurat_obj[['Mature_luminal']] <- PercentageFeatureSet(seurat_obj, features = f_lum)
  seurat_obj[['Luminal_progenitor']] <- PercentageFeatureSet(seurat_obj, features = f_prog)
  seurat_obj[['Luminal']] <- PercentageFeatureSet(seurat_obj, features = f_combined)
  seurat_obj[['percent.marker']] <- PercentageFeatureSet(seurat_obj, features = f_all)

  # Filter based on marker expression
  keep_cells <- WhichCells(seurat_obj, expression = Luminal > 0 | percent.marker > 0)
  seurat_obj <- subset(seurat_obj, cells = keep_cells)
  seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj), assay = "RNA")

  # Dimension reduction
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
  seurat_obj <- RunPCA(seurat_obj, features = f_combined)
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30)
  seurat_obj <- FindClusters(seurat_obj, resolution = 0.1)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)
  seurat_obj$orig.ident <- paste0("scCUTnTAG_", mark)

  assign(seurat_obj_name, seurat_obj)
}

#### Optional: Integration with scRNA ####
# scRNA <- readRDS("/path/to/scRNA_seurat_object.rds")
# scRNA_scCUTAG <- merge(scRNA, c(H3K4me3_seurat_filtered, H3K27ac_seurat_filtered))
# scRNA_scCUTAG <- NormalizeData(scRNA_scCUTAG, assay = "RNA")
# scRNA_scCUTAG <- ScaleData(scRNA_scCUTAG, features = features_all)
# scRNA_scCUTAG <- RunPCA(scRNA_scCUTAG, features = features_all)
# scRNA_scCUTAG <- FindNeighbors(scRNA_scCUTAG, dims = 1:30)
# scRNA_scCUTAG <- FindClusters(scRNA_scCUTAG, resolution = 0.1)
# scRNA_scCUTAG <- RunUMAP(scRNA_scCUTAG, dims = 1:10)

# anchors <- FindIntegrationAnchors(
#   object.list = list(scRNA, H3K4me3_seurat_filtered, H3K27ac_seurat_filtered),
#   anchor.features = features_all, reduction = "rpca", dims = 1:30, k.anchor = 100
# )

# scRNA_scCUTAG.integrated <- IntegrateEmbeddings(
#   anchorset = anchors,
#   reductions = scRNA_scCUTAG[["pca"]],
#   new.reduction.name = "integrated_pca",
#   dims.to.integrate = 1:30
# )

# scRNA_scCUTAG.integrated <- ScaleData(scRNA_scCUTAG.integrated)
# scRNA_scCUTAG.integrated <- FindNeighbors(scRNA_scCUTAG.integrated, reduction = "integrated_pca", dims = 1:15)
# scRNA_scCUTAG.integrated <- FindClusters(scRNA_scCUTAG.integrated, resolution = 0.1)
# scRNA_scCUTAG.integrated <- RunUMAP(scRNA_scCUTAG.integrated, reduction = "integrated_pca", dims = 1:10)

message("Gene activity, marker scoring, dimensional reduction, and integration setup completed.")