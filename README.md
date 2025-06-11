# IT-scC&T-seq Data Processing Pipeline

This repository contains a modular and reproducible pipeline for processing indexed tagmentation-based single-cell CUT&Tag-sequencing (IT-scCUT&Tag-seq) data. The workflow supports multiple histone marks and replicates, and integrates ArchR, Seurat + Signac for downstream single-cell analysis.

## 📁 Repository Structure
```
scCT-pipeline/
├── scripts/                # Shell-based preprocessing
│   ├── 01_demultiplex.sh  # Cutadapt-based index and linker trimming
│   ├── 02_alignment.sh    # Bowtie2 alignment
│   ├── 03_deduplication.sh # MarkDuplicates
│   ├── 04_filter_split.sh # Sinto-based barcode splitting
│   ├── 05_fragments.sh     # Generate fragment BED files
│   ├── 06_callpeaks.sh     # MACS2 peak calling
│
├── R/                     # R-based downstream analysis
│   ├── Signac-0.R         # Create and merge Seurat objects per histone mark for mammary gland histone mark datasets
│   ├── Signac-1.R         # Gene activity + marker scoring + UMAP for mammary gland histone mark datasets
│   ├── ArchR-0.R         # Create ArchR object and dimension reduction for single-cell LADs profiling
│
├── index/                 # Barcode FASTA files used by cutadapt
├── README.md              # This file
```

---

## 🛠️ Requirements
- bash ≥ 4.0
- `cutadapt`, `bowtie2`, `samtools`, `sambamba`, `macs2`, `picard`, `sinto`
- R ≥ 4.1.0 with the following packages:
  - `Seurat`, `Signac`, `Matrix`, `EnsDb.Mmusculus.v79`, `dplyr`, `GenomeInfoDb`, `GenomicRanges`, `ArchR`


---

## 📤 Citation
Please cite:
> Ma J#, Jin W#*, et al. IT-scCUT&Tag streamlines scalable, parallel profiling of protein-DNA interactions in single cells. *Genome Biology* (in revision).

---

For any questions or contributions, feel free to open an issue or pull request.

