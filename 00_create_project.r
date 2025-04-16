# 00_create_project.R

# User-defined parameters
PROJECT="ITscCT"
TISSUE="CL" #CellLines
MARK="LADs"
THREADS=10
# Paths and files
BASE_DIR <- "/path/to/project"
FRAG_DIR <- file.path(BASE_DIR, "fragments")
FRAG_FILE <- file.path(FRAG_DIR, paste0(PROJECT, "-", TISSUE, "-", Mark, ".sinto.sort.bed.gz"))

# Replicate annotation: sample sheet with barcode-to-replicate mapping
# Format: first column = cell barcode, second column = replicate name
REPLICATE_ANNOT <- file.path(BASE_DIR, "scLADs_workbook_replicates.csv")

# Optional: sample-to-celltype mapping (for future steps)
INDEX_MAP_FILE <- file.path(BASE_DIR, "index", paste0(PROJECT, "-", TISSUE, "-", TIMEPOINT, "-", REPLICATE, ".barcode_to_sample.txt"))

# ArchR setup
library(ArchR)
addArchRGenome("hg38")
addArchRThreads(threads = THREADS)
set.seed(26)

# Create output folder for initial project (proj-0)
OUTPUT_DIR_0 <- file.path(BASE_DIR, "ArchR", paste0(PROJECT, "-", TISSUE, "-", TIMEPOINT, "-", REPLICATE, "-proj-0"))
dir.create(OUTPUT_DIR_0, recursive = TRUE, showWarnings = FALSE)
setwd(OUTPUT_DIR_0)

# Input
inputFiles <- c("lads" = FRAG_FILE)

# Create Arrow files
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  filterTSS = 0,
  filterFrags = 0,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  TileMatParams = list(tileSize = 50000),
  minTSS = NULL,
  minFrags = NULL,
  threads = THREADS
)

# Create project (proj-0)
proj <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = OUTPUT_DIR_0,
  copyArrows = TRUE
)

# === Add cell type annotation from INDEX_MAP_FILE ===
barcode_map <- read.table(INDEX_MAP_FILE, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(barcode_map) <- c("CB", "CellType")
barcode_map$cellNames <- paste0("lad#", barcode_map$CB)

proj <- addCellColData(
  ArchRProj = proj,
  data = barcode_map$CellType,
  cells = barcode_map$cellNames,
  name = "cellType",
  force = TRUE
)

# === Add replicate information ===
replicate_map <- read.csv(REPLICATE_ANNOT, header = FALSE, stringsAsFactors = FALSE)
colnames(replicate_map) <- c("CB", "Replicate")
replicate_map$cellNames <- paste0("lad#", replicate_map$CB)

proj <- addCellColData(
  ArchRProj = proj,
  data = replicate_map$Replicate,
  cells = replicate_map$cellNames,
  name = "replicate",
  force = TRUE
)

# === Save enhanced project (proj-1) ===
OUTPUT_DIR_1 <- file.path(BASE_DIR, "ArchR", paste0(PROJECT, "-", TISSUE, "-", TIMEPOINT, "-", REPLICATE, "-proj-1"))
saveArchRProject(proj, outputDirectory = OUTPUT_DIR_1, overwrite = TRUE)

# 02_dimension_reduction.R

# User-defined parameters
PROJECT="ITscCT"
TISSUE="CL" #CellLines
MARK="LADs"
THREADS=10

# Directories
BASE_DIR <- "/path/to/project"
PROJ_INPUT_DIR <- file.path(BASE_DIR, "ArchR", paste0(PROJECT, "-", TISSUE, "-", MARK, "-proj-1"))
PROJ_OUTPUT_DIR <- file.path(BASE_DIR, "ArchR", paste0(PROJECT, "-", TISSUE, "-", MARK, "-proj-2"))

# Load ArchR
library(ArchR)
addArchRGenome("hg38")
addArchRThreads(THREADS)

# Load previous project (with metadata)
proj <- loadArchRProject(path = PROJ_INPUT_DIR)

# Run Iterative LSI on TileMatrix (default for LADs: 50kb bins)
proj <- addIterativeLSI(
  ArchRProj = proj,
  useMatrix = "TileMatrix",
  name = "IterativeLSI",
  iterations = 2,
  dimsToUse = 1:15,
  clusterParams = list(
    resolution = c(0.1, 0.2),
    sampleCells = 4000,
    n.start = 10
  ),
  outlierQuantiles = c(0.025, 0.975),
  verbose = TRUE,
  force = TRUE
)

# Add UMAP embedding
proj <- addUMAP(
  ArchRProj = proj,
  reducedDims = "IterativeLSI",
  name = "UMAP",
  nNeighbors = 20,
  minDist = 0.3,
  metric = "cosine",
  force = TRUE
)

# Add clustering using Seurat backend
proj <- addClusters(
  input = proj,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.2,
  force = TRUE,
  testBias = TRUE
)

# Save new project with LSI+UMAP+Clusters
saveArchRProject(
  ArchRProj = proj,
  outputDirectory = PROJ_OUTPUT_DIR,
  overwrite = TRUE
)