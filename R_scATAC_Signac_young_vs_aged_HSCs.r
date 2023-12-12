## load packages

library(Signac)
library(Seurat)
library(tidyverse)
library(patchwork) # plot | plot
library(GenomicRanges) # creating peak sets
library(future) # parallel processing
library(EnsDb.Mmusculus.v79) # annotation

# Suppress warnings
options(warn = -1)

# change the current plan to access parallelization
plan("multisession", workers = 28)
plan()

# increase size of plots from default
options(repr.plot.width = 14, 
        repr.plot.height = 14) # from 7, 7

# load peaks called on individual datasets
young_peaks <- read.table(
  file = "youngPeaks.txt",
  col.names = c("chr", "start", "end")
)
aged_peaks <- read.table(
  file = "agedPeaks.txt",
  col.names = c("chr", "start", "end")
)
young_peaks <- makeGRangesFromDataFrame(young_peaks)
aged_peaks <- makeGRangesFromDataFrame(aged_peaks)

# create common peak set
common_peaks <- reduce(x = c(young_peaks, aged_peaks))
common_peaks

# Filter out "bad" peaks based on width
peakwidths <- width(common_peaks)
common_peaks <- common_peaks[peakwidths  < 10000 & peakwidths > 20]
common_peaks

# load metadata
metadata_young <- read.table(
  file = "GSM5723631_Young_HSC_singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] # remove 1st row containing summary of frags not associated with any whitelisted barcodes
metadata_aged <- read.table(
  file = "GSM5723632_Aged_HSC_singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] # remove 1st containing summary of frags not associated with any whitelisted barcodes

head(metadata_young)
head(metadata_aged)

# perform an initial filtering of low count "cells" (barcodes)
nrow(metadata_young)
nrow(metadata_aged)

metadata_young <- metadata_young[metadata_young$passed_filters > 500, ]
metadata_aged <- metadata_aged[metadata_aged$passed_filters > 500, ]

nrow(metadata_young)
nrow(metadata_aged)

nrow(metadata_young)

# create fragment objects
frags_young <- CreateFragmentObject(
  path = "GSM5723631_Young_HSC_fragments.tsv.gz",
  cells = rownames(metadata_young)
)
frags_aged <- CreateFragmentObject(
  path = "GSM5723632_Aged_HSC_fragments.tsv.gz",
  cells = rownames(metadata_aged)
)
frags_young
frags_aged
str(frags_young)
str(frags_aged)

young_counts <- FeatureMatrix(
  fragments = frags_young,
  features = common_peaks,
  cells = rownames(metadata_young)
)

aged_counts <- FeatureMatrix(
  fragments = frags_aged,
  features = common_peaks,
  cells = rownames(metadata_aged)
) # ~6 minutes

str(young_counts)
dim(young_counts)
str(aged_counts)
dim(aged_counts)

young_assay <- CreateChromatinAssay(young_counts, fragments = frags_young)
young <- CreateSeuratObject(young_assay, assay = "ATAC", meta.data=metadata_young)

aged_assay <- CreateChromatinAssay(aged_counts, fragments = frags_aged)
aged <- CreateSeuratObject(aged_assay, assay = "ATAC", meta.data=metadata_aged)

young
aged

# add information to identify dataset of origin
young$dataset <- 'young'
aged$dataset <- 'aged'

# merge the datasets, adding a cell ID to make sure cell names are unique
combined <- merge(
  x = young,
  y = aged,
  add.cell.ids = c("young", "aged")
) # ~10 seconds

combined[["ATAC"]]

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)

# change to UCSC style 
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"

# add the gene information to the object
Annotation(combined) <- annotations

combined

# function to filter dataset on 5 QC metrics

qc_filter <- function(dataset){
    # compute QC metrics for each cell
    dataset <- NucleosomeSignal(object = dataset) # nucleosome signal score
    dataset <- TSSEnrichment(object = dataset, fast = FALSE) # TSS enrichment score
    dataset$pct_reads_in_peaks <- dataset$peak_region_fragments / dataset$passed_filters * 100 # fraction reads in peaks
    dataset$blacklist_ratio <- dataset$blacklist_region_fragments / dataset$peak_region_fragments # blacklist ratio

    # plot QC metrics PRE-filter
    print(
    VlnPlot(
      object = dataset,
      features = c('peak_region_fragments', 'TSS.enrichment', 'blacklist_ratio', 
                   'nucleosome_signal', 'pct_reads_in_peaks'),
      pt.size = 0.1,
      ncol = 5
    )) # print required to display
    
    # filter cells based on the 5 QC metrics (what cutoffs did the paper use? probs here are arbitrary)
    low_prf <- quantile(dataset[["peak_region_fragments"]]$peak_region_fragments, probs = 0.025)
    high_prf <- quantile(dataset[["peak_region_fragments"]]$peak_region_fragments, probs = 0.975)
    low_prp <- quantile(dataset[["pct_reads_in_peaks"]]$pct_reads_in_peaks, probs = 0.05)
    high_blr <- quantile(dataset[["blacklist_ratio"]]$blacklist_ratio, probs = 0.95)
    high_ns <- quantile(dataset[["nucleosome_signal"]]$nucleosome_signal, probs = 0.95)
    low_ts <- quantile(dataset[["TSS.enrichment"]]$TSS.enrichment, probs = 0.05)

    dataset <- subset(
                  x = dataset,
                  subset = peak_region_fragments > low_prf &
                            peak_region_fragments < high_prf &
                            pct_reads_in_peaks > low_prp &
                            blacklist_ratio < high_blr &
                            nucleosome_signal < high_ns &
                            TSS.enrichment > low_ts
    )
    return(dataset)    
}

# do QC filter to exclude low quality cells
combined
combined <- qc_filter(combined) # ~200s
combined

# plot QC metrics POST-filter
VlnPlot(
  object = combined,
  features = c('peak_region_fragments', 'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal', 'pct_reads_in_peaks'),
  pt.size = 0.1,
  ncol = 5
)

# Normalization (TF-IDF) and linear dimensional reduction (SVD); TF-IDF + SVD = LSI
combined <- RunTFIDF(combined) %>% 
              FindTopFeatures(min.cutoff = 20) %>% 
                RunSVD()
combined

# Graph-based clustering & non-linear dimension reduction for visualization
combined <- RunUMAP(object = combined, reduction = 'lsi', dims = 2:30) %>%
              FindNeighbors(reduction = 'lsi', dims = 2:30) %>% 
                FindClusters(verbose = TRUE, algorithm = 3)


# plot UMAP
DimPlot(object = combined, 
        pt.size = 3, label.size = 10, 
        shuffle = TRUE, 
        label = TRUE) + NoLegend()

# UMAP colored by young vs aged
DimPlot(object = combined, 
        pt.size = 3, label.size = 10, 
        shuffle = TRUE, 
        label = TRUE, 
       group.by = "dataset") + NoLegend()

# aged vs young per cluster
table(combined$dataset, combined$ATAC_snn_res.0.8)

gene.activities <- GeneActivity(combined) # ~ 130s

gene.activities[1:10, 1:66]

# Explicitly close multisession workers by switching plan - otherwise NormalizeData below crashes
plan(sequential)
plan()

# add the gene activity matrix to the Seurat object as a new assay and normalize it
combined[['RNA']] <- CreateAssayObject(counts = gene.activities)

combined <- NormalizeData(
  object = combined,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(combined$nCount_RNA)
) # ~3s

combined[['RNA']]

# visualize the activities of canonical marker genes
DefaultAssay(combined) <- 'RNA'

plot1 <- FeaturePlot(
          object = combined,
          features = c('Cd34', 'Cd38', 'Thy1', 'Ptprc', 
                       'Prom1', 'Kit', 'Ly6a', 'Flt3', 
                       'Slamf1', 'Pecam1', 'Itgam'),
          pt.size = 0.1,
          max.cutoff = 'q95',
          ncol = 3
        )
plot1

# change back to working with peaks instead of gene activities
Idents(combined) <- 'dataset'
DefaultAssay(combined) <- 'ATAC'

# find differential markers for the two identity classes (young vs aged)
da_peaks <- FindMarkers(
  object = combined,
  ident.1 = "young",
  ident.2 = "aged",
  test.use = 'LR', # use logistic regression for differential accessibility test
    min.pct = 0.05, # only consider peaks found in at least 5% of the cells of either of the two groups
  latent.vars = 'peak_region_fragments' # latent variable to mitigate effect of differential sequencing depth on result 
)

head(da_peaks)

length(da_peaks)

length(da_peaks[da_peaks$p_val_adj < 0.001, ])

sum(da_peaks[da_peaks$p_val_adj < 0.001, ]$avg_log2FC > 0)

# visualization
plot1 <- VlnPlot(
  object = combined,
  features = rownames(da_peaks)[1],
  pt.size = 2,
  idents = c("young","aged")
)
plot2 <- FeaturePlot(
  object = combined,
  features = rownames(da_peaks)[1],
  pt.size = 2
)

plot1 | plot2

rownames(da_peaks)[1]

# Sort the data frame by the absolute value of avg_log2FC
da_peaks_sorted <- da_peaks[order(-abs(da_peaks$avg_log2FC)), ]
head(da_peaks_sorted)

da_peaks_sorted[da_peaks_sorted$avg_log2FC > 0, ] %>% head

# visualization
plot1 <- VlnPlot(
  object = combined,
  features = 'chr19-6290930-6294154',
  pt.size = 2,
  idents = c("young","aged")
)
plot2 <- FeaturePlot(
  object = combined,
  features = 'chr19-6290930-6294154',
  pt.size = 2
)

plot1 | plot2

da_peaks_closest_feature <- ClosestFeature(combined, 
                                          regions = rownames(da_peaks))
head(da_peaks_closest_feature)

da_peaks_closest_feature[da_peaks_closest_feature$distance == 0, ] %>% head(n = 25)

CoveragePlot(
  object = combined,
  region = rownames(da_peaks)[1],
  extend.upstream = 1000,
  extend.downstream = 5000
)

CoveragePlot(
  object = combined,
  region = "chr19-20599880-20607776",
  extend.upstream = 1000,
  extend.downstream = 5000
)

CoveragePlot(
  object = combined,
  region = "chr5-26903642-26905903",
  extend.upstream = 1000,
  extend.downstream = 5000
)

CoveragePlot(
  object = combined,
  region = "chr1-15285819-15288023",
  extend.upstream = 1000,
  extend.downstream = 5000
)

sessionInfo()

