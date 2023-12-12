# scATAC-seq: Comparing chromatin accessibility in young vs aged hematopoietic stem cells (HSCs) in young vs aged mice

I wrote this R/Bioconductor script analysis pipeline for single-cell ATAC-seq (scATAC-seq) data looking at differential accessibility in hematopoietic stem cells (HSCs) from 10-week-old 'young' mice vs 20-month-old 'aged' mice. [The data](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE162662) is from the NCBI Gene Expression Omnibus (GEO) repository.

### Overall Functionality:
1. **Data Loading and Preprocessing**: It loads necessary R packages (Signac, Seurat, tidyverse, etc.) and reads in peak files, metadata, and fragment files for young and aged single-cell ATAC-seq datasets.
2. **Peak Filtering and Common Set Creation**: It identifies a common set of peaks by reducing peaks from individual datasets and filters out peaks based on width criteria.
3. **Initial Quality Control**: Filters out low-quality cells based on specific cutoffs for various quality control metrics.
4. **Count Matrix Generation**: Generates count matrices for both young and aged datasets based on the common peak set.
5. **Seurat Object Creation**: Constructs Seurat objects for each dataset.
6. **Integration and Dimensionality Reduction**: Integrates datasets, performs TF-IDF normalization followed by SVD, and conducts UMAP-based visualization and clustering.
7. **Gene Annotation and Analysis**: Extracts gene annotations, adds them to the Seurat object, and conducts gene activity analysis.
8. **Normalization and Visualization of "RNA" Data**: Normalizes gene activity "RNA" data, visualizes canonical marker genes, and identifies differential peaks between young and aged datasets.
9. **Visualization of Peaks**: Visualizes coverage plots for selected peaks and their closest genomic features.

### Input Files Required:
- Individual peak files for young and aged datasets (`youngPeaks.txt`, `agedPeaks.txt`)
- Metadata files for young and aged datasets (`GSM5723631_Young_HSC_singlecell.csv`, `GSM5723632_Aged_HSC_singlecell.csv`)
- Fragment files for young and aged datasets (`GSM5723631_Young_HSC_fragments.tsv.gz`, `GSM5723632_Aged_HSC_fragments.tsv.gz`)
- Annotation file (`EnsDb.Mmusculus.v79`)

### R/Bioconductor Packages and Bioinformatics Tools Required:
- Signac
- Seurat
- tidyverse
- patchwork
- GenomicRanges
- future
- EnsDb.Mmusculus.v79 (annotation)


### Outputs:
- Seurat objects (`combined`, `young`, `aged`) containing integrated and processed single-cell ATAC-seq data.
- Visualizations: Various plots for quality control, gene activity, differential peaks, UMAP visualization, coverage plots for peaks, and more.
