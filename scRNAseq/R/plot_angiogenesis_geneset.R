suppressMessages(source("R/load_data_scran.R"))

### functions

trajectory.inference <- function(x, ... ){
  require(slingshot)
  require(mclust)
  require(magrittr)
  x %<>% runDiffusionMap(...)
  cl <- Mclust(reducedDim(x, "DiffusionMap"))$classification
  x %<>% slingshot(clusterLabels = cl, reducedDim = "DiffusionMap")
  x$slingPseudotime_all <- rowMeans(as.data.frame(colData(x)[,grepl("slingPseudotime",colnames(colData(x)))]), na.rm = T)
  x
}

### end

### infer pseudotime trajectory
sce.nk.ss2.filt %<>% trajectory.inference(feature_set = head(hvgs.ss2,500))

## load angiogenesis GO term genes from http://amigo.geneontology.org/amigo/term/GO:0001525 (accessed: 2020-11-01)
go.angiogenesis <- fread("GO/GO0001525_angiogenesis.txt",header = F)

## subset to angiogenesis genes
dat.angio <- assay(sce.nk.ss2.filt, "scaled_logcounts")[gsub("_.*","",rownames(sce.nk.ss2.filt)) %in% go.angiogenesis$V1,]
# keep genes detected in >10% of cells
dat.angio.filt <- dat.angio[which(rowMeans(dat.angio > 0,na.rm = T) > 0.1),]

## melt data and add annotations
dat.angio.melt <- as.data.table(melt(t(dat.angio.filt)))
#pseudotime
dat.angio.melt[, pseudotime := rep(sce.nk.ss2.filt$slingPseudotime_all, length(unique(Var2)))]
# CD39/KIR status
# find KIRs (not KIR2DL4)
idx.kir <- grep("KIR[[:alnum:]]D.[1,2,3]",rownames(sce.nk.ss2.filt))
dat.angio.melt[,group := rep(paste0("KIR", c("-","+")[(colSums(logcounts(sce.nk.ss2.filt)[idx.kir,])>0) + 1], "CD39",c("-","+")[(logcounts(sce.nk.ss2.filt)["ENTPD1_ENSG00000138185",]>0) + 1]), length(unique(Var2)))]

# calculate pseudotime correlation
dat.angio.cor <- dat.angio.melt[, mean(value, na.rm=T), by=pseudotime][,cor.test(pseudotime, V1, method="spearman", exact = F)]

## export correlations
# fwrite(dat.angio.melt[,list(r = cor(value,pseudotime), rho = cor(value,pseudotime, method="spearman") ), by="Var2"][order(-r,-rho), list(gene=Var2, r, rho)], "data/angiogenesis_genes_correlations.tsv", quote = F, sep = "\t")

# p-values
dat.angio.melt[!group == "KIR-CD39+", mean(value, na.rm=T), by=c("Var1","group")][,pairwise.wilcox.test(V1, group, p.adjust.method = "fdr")]
dat.angio.melt[!group == "KIR-CD39+", mean(value, na.rm=T), by=c("Var1","group")][,TukeyHSD(aov(V1~group))]$group

## plot
library(ggplot2)
library(cowplot)
# pseudotime
p.angio <-
ggplot(dat.angio.melt[, mean(value, na.rm=T), by=pseudotime], aes(x=pseudotime, y=V1, col=pseudotime)) +
  geom_density2d(col="black", alpha=0.5) +
  geom_point() +
  annotate("text", x = 0, y = 0.35, hjust=0, label=paste0("rho = ", round(dat.angio.cor$estimate,3),"\n","p = ",signif(dat.angio.cor$p.value,3))) +
  labs(x="Pseudotime", y="Average Z-score", title="GO:0001525 Angiogenesis") +
  coord_cartesian(ylim=c(-0.4,0.4)) +
  expand_limits(y=0) +
  scale_colour_viridis_c(option = "plasma") +
  theme_cowplot()

ggsave2("plots/scatter_pseudotime_v_angiogenesis.pdf", width = 6, height = 4, p.angio)

# group
p.angio.box <-
ggplot(dat.angio.melt[!group == "KIR-CD39+", mean(value, na.rm=T), by=c("Var1","group")], aes(x=group, y=V1, col=group)) +
  geom_violin(scale="width", show.legend = F) +
  ggbeeswarm::geom_quasirandom(show.legend = F) +
  geom_hline(yintercept = 0, lty=2) +
  labs(x=NULL, y="Average Z-score", title="GO:0001525 Angiogenesis") +
  coord_cartesian(ylim=c(-0.4,0.4)) +
  scale_color_manual(values = rev(pal.qual) ) +
  expand_limits(y=0) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))

ggsave2("plots/violin_group_v_angiogenesis_quasirandom.pdf", width = 4, height = 4, p.angio.box)
