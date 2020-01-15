suppressMessages(source("R/load_data.R"))

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

geom_slingshot <- function(data, ... ){
  require(slingshot)
  require(ggplot2)
  require(data.table)
  scs <- slingCurves(data)
  pc <- rbindlist(lapply(slingCurves(data),function(x) with(x,data.table(s[ord,1:2])) ),idcol = T)
  return(geom_path(data = pc, aes(DC1, DC2, group=.id), ... ))
}

reformat.zlm <- function(fit, hypothesis){
  require(MAST)
  res <- summary(fit, doLRT=hypothesis)
  res.dt <- res$datatable
  hurdle <- merge(res.dt[contrast==hypothesis & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                  res.dt[contrast==hypothesis & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
  hurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
  setorder(hurdle,fdr)
  return(hurdle)
}

plotHeat <- function(x, features, exprs_values = "scaled_logcounts", rowlabels, cutoff = 0.1, ... ){
  require(scater)
  require(ComplexHeatmap)
  require(circlize)
  require(RColorBrewer)
  require(viridis)
  
  dat <- na.omit(assay(x,exprs_values)[features,])
  dat <- subset(dat, rowMeans(dat>0)>cutoff)
  rha <- rowAnnotation(foo = anno_mark(at = match(rowlabels,row.names(dat)), labels = gsub("_.*","",rowlabels)))
  cha <- with(colData(x),HeatmapAnnotation(Pseudotime=slingPseudotime_all, col=list(Pseudotime=colorRamp2(seq(0,0.4,length.out = 5),viridis_pal(option="plasma")(5))),annotation_legend_param = list(Pseudotime = list(direction = "horizontal"))))
  
  Heatmap(dat,
    show_column_names = F,
    show_row_names = F,
    show_row_dend = F,
    cluster_columns = F,
    column_order = order(x$slingPseudotime_all),
    name = "Z-scores",
    col= colorRamp2(seq(-2,2,length.out = 11),rev(brewer.pal(11,"RdBu"))),
    border="black",
    top_annotation = cha,
    right_annotation = rha,
    use_raster = T,
    raster_device = "CairoPNG",
    heatmap_legend_param = list(direction="horizontal"),
    ...
    )
}

### end

### infer pseudotime trajectory
sce.nk.ss2.filt %<>% trajectory.inference(feature_set = head(hvgs.ss2,500))

### plot diffusion maps
## plot pseudotime
# NOTE: diffusion maps can have a non-deterministic behaviour and plots can appear mirrored, see https://github.com/theislab/destiny/issues/5
plotDiffusionMap(sce.nk.ss2.filt, colour_by="slingPseudotime_all") + geom_slingshot(sce.nk.ss2.filt) + labs(x="DC1",y="DC2") + theme(legend.position = "top") + scale_fill_viridis_c(name="Pseudotime", option="inferno")

## plot gene expression over diffusion map
gnls <- c("CXCR4_ENSG00000121966", "CD44_ENSG00000026508","GZMA_ENSG00000145649", "ENTPD1_ENSG00000138185","KIR2DL1_ENSG00000125498", "LGALS1_ENSG00000100097","LGALS9_ENSG00000168961")

multiplot(cols = 4,
  plotlist = lapply(gnls, function(x) plotDiffusionMap(sce.nk.ss2.filt, colour_by=x) + geom_slingshot(sce.nk.ss2.filt) + labs(x="DC1",y="DC2") + theme(legend.position = "top") + scale_fill_viridis_c(name=gsub("_.*","",x)) )
)

### plot heatmap
## combined with surface proteome- and RNA-seq results
# load surface proteome and RNA-seq DEG gene lists
gnls.surf <- rownames(sce.nk.ss2.filt)[gsub("_.*","",rownames(sce.nk.ss2.filt)) %in% unlist(read.delim("genelists/genes_surface_proteome.txt",header = F))]
gnls.deg <- na.omit(rownames(sce.nk.ss2.filt)[match(unlist(read.delim("genelists/genes_deg.txt",header = F)), gsub(".*_","",rownames(sce.nk.ss2.filt)))])
gnls.degsurf <- unique(c(gnls.deg,gnls.surf))

# set up plot annotations and draw plot
rha.degsurf <- rowAnnotation(Source=cbind("Bulk RNAseq"=gnls.degsurf %in% gnls.deg, "Surface proteome" = gnls.degsurf %in% gnls.surf), border=T, annotation_legend_param = list(border="black"),col = list(Source=c("TRUE"="#0B0B0B","FALSE"="#F7F7F7")))
p.heat.degsurf.split <- plotHeat(sce.nk.ss2.filt, features = gnls.degsurf, rowlabels = head(gnls.degsurf[order(rowVars(logcounts(sce.nk.ss2.filt)[gnls.degsurf,]),decreasing = T)],50), row_split=2, row_title=NULL)

draw(rha.degsurf+p.heat.degsurf.split,heatmap_legend_side="top",annotation_legend_side="top",merge_legend=T)
