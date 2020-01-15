library(data.table)
library(magrittr)

# colour palettes
pal.qual <- c("#FF9E4A", "#AD8BC9", "#A2A2A2")

## SS2
# load pre-processed data
meta.ss2 <- fread("E-MTAB-6678/meta_ss2.txt")
colnames(meta.ss2)[1] <- "sample"
dat.ss2 <- fread("E-MTAB-6678/raw_data_ss2.txt")
idx.nk.ss2 <- grepl("^dNK[[:alnum:]]",meta.ss2$annotation)

# load in scran
library(scater)
library(scran)
sce.nk.ss2 <- SingleCellExperiment(
  assays = list(
    counts = as.matrix(data.frame(dat.ss2[,c("Gene",meta.ss2[idx.nk.ss2,sample]), with = F],row.names = 1))
  ),
  colData = meta.ss2[idx.nk.ss2]
)

# remove outliers based on features detected and % in top features
sce.nk.ss2 %<>% calculateQCMetrics()
drop.features <- sce.nk.ss2$total_features_by_counts < 2e3
drop.pcttop <- sce.nk.ss2$pct_counts_in_top_500_features > 75
sce.nk.ss2$outlier <- drop.features | drop.pcttop
sce.nk.ss2.filt <- subset(sce.nk.ss2,,!outlier)

# cluster and normalise per bin
sce.nk.ss2.filt$cluster <- quickCluster(sce.nk.ss2.filt,use.ranks=F)
sce.nk.ss2.filt %<>% computeSumFactors(clusters=sce.nk.ss2$cluster)
sce.nk.ss2.filt %<>% normalize()

# identify hvgs
var.fit.ss2 <- trendVar(sce.nk.ss2.filt, use.spikes=F)
var.decomp.ss2 <- decomposeVar(sce.nk.ss2.filt, var.fit.ss2)
var.decomp.ss2$name <- rownames(var.decomp.ss2)
var.decomp.ss2 <- var.decomp.ss2[with(var.decomp.ss2,order(-bio,FDR)),]
hvgs.ss2 <- var.decomp.ss2[var.decomp.ss2$FDR < 0.05,'name']

# add scaled counts (z-score) to assays
assays(sce.nk.ss2.filt)[["scaled_logcounts"]] <- t(apply(assays(sce.nk.ss2.filt)[["logcounts"]],1,scale))
