library(SnapATAC);
library(GenomicRanges);
library(ggplot2);

##For Multiple Samples
snap.files = system("ls *.snap", intern = T)

sample.names = sapply(strsplit(snap.files, '\\.'), '[', 1)

barcode.files = system("ls ./singlecell_*.csv", intern = T);

x.sp.ls = lapply(seq(snap.files), function(i){
  createSnap(
    file=snap.files[i],
    sample=sample.names[i]
  );
})
names(x.sp.ls) = sample.names;

barcode.ls = lapply(seq(snap.files), function(i){
  barcodes = fread(
    barcode.files[i], 
    header=TRUE
  );
  # remove NO BAROCDE line
  barcodes = barcodes[2:nrow(barcodes),];
  barcodes$logUMI = log10(barcodes$passed_filters + 1);
  barcodes$promoter_ratio = (barcodes$promoter_region_fragments+1) / (barcodes$passed_filters + 1);
  barcodes
})

plots = lapply(seq(snap.files), function(i){
  p1 = ggplot(
    barcode.ls[[i]], 
    aes(x=logUMI, y=promoter_ratio)) + 
    geom_point(size=0.3, col="grey") +
    theme_classic()	+
    ggtitle(sample.names[[i]]) +
    ylim(0, 1) + xlim(0, 6) + 
    labs(x = "log10(UMI)", y="promoter ratio")
  p1
})

pdf("covplots.pdf")
plots
dev.off()

cutoff.logUMI.low = array(3, length(barcode.ls));
cutoff.logUMI.high = array(5, length(barcode.ls));
cutoff.FRIP.low = array(0.125, length(barcode.ls));
cutoff.FRIP.high = array(0.5, length(barcode.ls));
barcode.ls = lapply(seq(snap.files), function(i){
  barcodes = barcode.ls[[i]];
  idx = which(
    barcodes$logUMI >= cutoff.logUMI.low[i] & 
      barcodes$logUMI <= cutoff.logUMI.high[i] & 
      barcodes$promoter_ratio >= cutoff.FRIP.low[i] &
      barcodes$promoter_ratio <= cutoff.FRIP.high[i]
  );
  barcodes[idx,]
});

x.sp.ls = lapply(seq(snap.files), function(i){
  barcodes = barcode.ls[[i]];
  x.sp = x.sp.ls[[i]];
  barcode.shared = intersect(x.sp@barcode, barcodes$barcode);
  x.sp = x.sp[match(barcode.shared, x.sp@barcode),];
  barcodes = barcodes[match(barcode.shared, barcodes$barcode),];
  x.sp@metaData = barcodes;
  x.sp
})
names(x.sp.ls) = sample.names

# combine snap objects
x.sp = Reduce(snapRbind, x.sp.ls);
x.sp@metaData["sample"] = x.sp@sample;
x.sp
table(x.sp@sample)

x.sp = addBmatToSnap(x.sp, bin.size=5000)
x.sp = makeBinary(x.sp, mat="bmat")

black_list = read.table("/media/RND/HDD-7/RZiffra/rziffra/hg38/hg38.blacklist.bed.gz");
black_list.gr = GRanges(
  black_list[,1], 
  IRanges(black_list[,2], black_list[,3])
);
idy = queryHits(
  findOverlaps(x.sp@feature, black_list.gr)
);
if(length(idy) > 0){
  x.sp = x.sp[,-idy, mat="bmat"];
};
x.sp

chr.exclude = seqlevels(x.sp@feature)[grep("random|chrM", seqlevels(x.sp@feature))];
idy = grep(paste(chr.exclude, collapse="|"), x.sp@feature);
if(length(idy) > 0){
  x.sp = x.sp[,-idy, mat="bmat"]
};
x.sp

bin.cov = log10(Matrix::colSums(x.sp@bmat)+1);
hist(
  bin.cov[bin.cov > 0], 
  xlab="log10(bin cov)", 
  main="log10(Bin Cov)", 
  col="lightblue", 
  xlim=c(0, 5)
);
bin.cutoff = quantile(bin.cov[bin.cov > 0], 0.95);
idy = which(bin.cov <= bin.cutoff & bin.cov > 0);
x.sp = x.sp[, idy, mat="bmat"];
x.sp

idx = which(Matrix::rowSums(x.sp@bmat) > 500);
x.sp = x.sp[idx,];
x.sp

x.sp = runLSA(
  obj = x.sp,
  input.mat = "bmat",
  pc.num = 50,
  logTF = TRUE
)

plotDimReductPW(
  obj=x.sp, 
  eigs.dims=1:50,
  point.size=0.3,
  point.color="grey",
  point.shape=19,
  point.alpha=0.6,
  down.sample=5000,
  pdf.file.name=NULL, 
  pdf.height=7, 
  pdf.width=7
);

x.sp = runKNN(
  obj=x.sp,
  eigs.dims=1:20,
  k=15
);

library(leiden)
x.sp=runCluster(
  obj=x.sp,
  tmp.folder=tempdir(),
  louvain.lib="leiden",
  resolution = 0.6,
  seed.use=10
);
x.sp@metaData$cluster = x.sp@cluster;

x.sp = runViz(
  obj=x.sp, 
  tmp.folder=tempdir(),
  dims=2,
  eigs.dims=1:20, 
  method="umap",
  seed.use=10
);

area_vector <-sapply(strsplit(x.sp@sample,"_"), `[`, 1)
specimen_vector <-sapply(strsplit(x.sp@sample,"_"), `[`, 2)

##Run scAlign batch correction. See sca.R script

par(mfrow = c(1, 1));
png(filename = "Fig2a_ATAConly_Final.png", height = 3000, width = 3000, res = 300)
plotViz(
  obj=x.sp,
  method="umap", 
  main="V1 ATAC Clusters",
  cex.main=2,
  point.color=x.sp@metaData$cluster, 
  point.size=0.5, 
  point.shape=20, 
  point.alpha=0.8, 
  text.add=FALSE,
  text.size=2,
  text.color="black",
  text.halo.add=TRUE,
  text.halo.color="white",
  text.halo.width=0.2,
  down.sample=80000,
  legend.add=F,
  legend.text.size = 1,
  legend.pos = "topright"
);

png(filename = "AllPrimary_Areas_Final.png", height = 3000, width = 3000, res = 300)
plotViz(
  obj=x.sp,
  method="umap", 
  main="All Primary Areas",
  point.color=area_vector, 
  point.size=0.2, 
  point.shape=20, 
  point.alpha=0.8, 
  text.add=F,
  text.size=0.75,
  text.color="black",
  text.halo.add=TRUE,
  text.halo.color="white",
  text.halo.width=0.1,
  down.sample=80000,
  legend.add=T,
  legend.text.size = 1,
  legend.pos = "topright"
);
png(filename = "AllPrimary_ReadDepth.png", height = 960, width = 960)
plotFeatureSingle(
  obj=x.sp,
  feature.value=log(x.sp@metaData[,"passed_filters"]+1,10),
  method="umap", 
  main="Log Read Depth",
  point.size=0.2, 
  point.shape=19, 
  down.sample=80000,
  quantiles=c(0.01, 0.99)
); 
png(filename = "AllPrimary_FRiP.png", height = 960, width = 960)
plotFeatureSingle(
  obj=x.sp,
  feature.value=x.sp@metaData$peak_region_fragments / x.sp@metaData$passed_filters,
  method="umap", 
  main="Fraction of Reads in Peaks",
  point.size=0.2, 
  point.shape=19, 
  down.sample=80000,
  quantiles=c(0.01, 0.99) # remove outliers
);
plotFeatureSingle(
  obj=x.sp,
  feature.value=x.sp@metaData$duplicate / x.sp@metaData$total,
  method="umap", 
  main="10X Brain Duplicate",
  point.size=0.2, 
  point.shape=19, 
  down.sample=10000,
  quantiles=c(0.01, 0.99) # remove outliers
);

genes = read.table("/media/RND/HDD-7/RZiffra/rziffra/hg38/hg38_gencode_genes.bed");


library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)

gene.coords <- genes(EnsDb.Hsapiens.v86, filter = ~ gene_biotype == "protein_coding")
seqlevelsStyle(gene.coords) <- 'UCSC'
genebody.coords <- keepStandardChromosomes(gene.coords, pruning.mode = 'coarse')
genebodyandpromoter.coords <- Extend(x = genebody.coords, upstream = 2000, downstream = 0)
genebodypromoter.df = as.data.frame(genebodyandpromoter.coords)
genebodypromoter.df = genebodypromoter.df[,c(1:3,7)]
genebodypromoter.df = genebodypromoter.df[grep("chrM",genebodypromoter.df[,1], invert = T),]

genes = genebodypromoter.df

genes.gr = GRanges(genes[,1], 
                     IRanges(genes[,2], genes[,3]), name=genes[,4]
);

# re-add the cell-by-bin matrix to the snap object;
x.sp = addBmatToSnap(x.sp);
x.sp = createGmatFromMat(
  obj=x.sp, 
  input.mat="bmat",
  genes=genes.gr,
  do.par=TRUE,
  num.cores=16
);
# normalize the cell-by-gene matrix
x.sp = scaleCountMatrix(
  obj=x.sp, 
  cov=x.sp@metaData$passed_filters + 1,
  mat="gmat",
  method = "log"
);
# smooth the cell-by-gene matrix
myRunMagic <- function (obj, input.mat, step.size) {
  A = obj@graph@mat;
  data.use = obj@gmat;
  
  # smooth
  A = A + t(A);
  A = A / Matrix::rowSums(A);
  data.use.smooth = A %*% data.use;
  if(step.size > 1){
    for(i in 1:step.size){
      data.use.smooth = A %*% data.use.smooth;
    }
  }
  
  slot(obj, input.mat) = data.use.smooth;    
  return(obj)
}

x.sp = myRunMagic(
  obj=x.sp,
  input.mat="gmat",
  step.size=3
);

marker.genes = c(
  "HES1","SOX2","GFAP","SATB2","NEUROD2","EOMES","FEZF2","TBR1",
  "BHLHE22","HOPX","CRYAB","AIF1","OLIG1","AQP4","DLX1","LHX6")

genes.sel.gr <- genes.gr[which(genes.gr$name %in% marker.genes)];

par(mfrow = c(4, 4));
for(i in 1:length(marker.genes)){
  png(paste0("AllOrg_",marker.genes[i],".png"), height = 3000, width = 3000, res = 300)
  plotFeatureSingle(
    obj=x.sp,
    feature.value=x.sp@gmat[, marker.genes[i]],
    method="umap", 
    main=marker.genes[i],
    point.size=0.2, 
    point.shape=19, 
    down.sample=50000,
    quantiles=c(0.01,0.99)
  )
  dev.off()  
};

ensemble.ls = lapply(split(seq(length(x.sp@cluster)), x.sp@cluster), function(x){
  SnapATAC::colMeans(x.sp[x,], mat="bmat");
})
# cluster using 1-cor as distance  
hc = hclust(as.dist(1 - cor(t(do.call(rbind, ensemble.ls)))), method="ward.D2");
par(mfrow=c(1,2))
plotViz(
  obj=x.sp,
  method="umap", 
  main="10X Brain Cluster",
  point.color=x.sp@cluster, 
  point.size=0.2, 
  point.shape=19, 
  point.alpha=0.8, 
  text.add=TRUE,
  text.size=0.75,
  text.color="black",
  text.halo.add=TRUE,
  text.halo.color="white",
  text.halo.width=0.1,
  down.sample=30000,
  legend.add=TRUE,
  legend.pos = "right",
  legend.text.size = 0.5
);
plot(hc, hang=-1, xlab="");


peaks.gr = runMACSForAll(
  obj=x.sp,
  path.to.snaptools="/home/rziffra/.local/bin/snaptools",
  path.to.macs="/usr/local/bin/macs2",
  output.prefix="AllPrimary_downsampled",
  num.cores=16,
  min.cells=50,
  gsize="hs", 
  buffer.size=500, 
  macs.options="--nomodel --shift -37 --ext 73 --qval 5e-2 -B --SPMR",
  tmp.folder="/media/RND/HDD-7/RZiffra/rziffra/allsnap/primary/temp/"
)

x.sp = createPmat(
  x.sp, 
  peaks=peaks.gr,
  ncell.chunk=20,
  do.par=TRUE,
  num.cores=24
)

overlap_peaks = list()
for(i in 1:length(peak.gr.ls)){overlap_peaks[[i]] = peak.gr[unique(findOverlaps(peak.gr.ls[[i]],peak.gr, select = "first")),]}

peaks.df = as.data.frame(peak.gr)[,1:3];
write.table(peaks.df,file = "peaks.gw17.bed",append=FALSE,
              quote= FALSE,sep="\t", eol = "\n", na = "NA", dec = ".", 
              row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"),
              fileEncoding = "")
saveRDS(x.sp, file="gw17.snap.rds");

x.sp = readRDS("gw17.snap.rds");
x.sp = addPmatToSnap(x.sp);
x.sp = makeBinary(x.sp, mat="pmat");
x.sp

DARs = findDAR(
  obj=x.sp,
  input.mat="pmat",
  cluster.pos=7,
  cluster.neg = NULL,
  cluster.neg.method="knn",
  test.method="exactTest",
  bcv=0.2, #0.4 for human, 0.1 for mouse
  seed.use=10
);
DARs$FDR = p.adjust(DARs$PValue, method="BH");
idy = which(DARs$FDR < 5e-2 & DARs$logFC > 0);
head(DARs[order(DARs$FDR),])
par(mfrow = c(1, 2));
plot(DARs$logCPM, DARs$logFC, 
       pch=19, cex=0.1, col="grey", 
       ylab="logFC", xlab="logCPM",
       main="Top 1000 Peaks - Cluster 7"
)
points(DARs$logCPM[idy], 
         DARs$logFC[idy], 
         pch=19, 
         cex=0.5, 
         col="red"
)
abline(h = 0, lwd=1, lty=2);
covs = Matrix::rowSums(x.sp@pmat);
vals = Matrix::rowSums(x.sp@pmat[,idy]) / covs;
vals.zscore = (vals - mean(vals)) / sd(vals);
png("tRG_peaks.png",res = 300,height = 3000, width = 3000)
plotFeatureSingle(
  obj=x.sp,
  feature.value=vals.zscore,
  method="umap", 
  main="tRG Specific Peaks",
  point.size=0.2, 
  point.shape=19, 
  down.sample=25000,
  quantiles=c(0.01, 0.99)
)
dev.off()
idy.ls = lapply(levels(x.sp@cluster), function(cluster_i){
  DARs = findDAR(
    obj=x.sp,
    input.mat="pmat",
    cluster.pos=cluster_i,
    cluster.neg=NULL,
    cluster.neg.method="knn",
    bcv=0.2,
    test.method="exactTest",
    seed.use=10
  );
  DARs$FDR = p.adjust(DARs$PValue, method="BH");
  idy = which(DARs$FDR < 5e-2 & DARs$logFC > 0);
  if((x=length(idy)) < 1000L){
    PValues = DARs$PValue;
    PValues[DARs$logFC < 0] = 1;
    idy = order(PValues, decreasing=FALSE)[1:1000];
    rm(PValues); # free memory
  }
  idy
})
names(idy.ls) = levels(x.sp@cluster);
par(mfrow = c(1,1));
for(cluster_i in levels(x.sp@cluster)){
  print(cluster_i)
  idy = idy.ls[[cluster_i]];
  vals = Matrix::rowSums(x.sp@pmat[,idy]) / covs;
  vals.zscore = (vals - mean(vals)) / sd(vals);
  plotFeatureSingle(
    obj=x.sp,
    feature.value=vals.zscore,
    method="umap", 
    main=cluster_i,
    point.size=0.1, 
    point.shape=19, 
    down.sample=5000,
    quantiles=c(0.01, 0.99)
  );
}

clusterpeaks = list()
for(i in 1:length(idy.ls)){clusterpeaks[[i]]=peak.gr[idy.ls[[i]],]}

motifs = list()
for(i in 1:length(idy.ls)){
  motifs[[i]] = runHomer(
  x.sp[,idy.ls[[i]],"pmat"], 
  mat = "pmat",
  path.to.homer = "/home/rziffra/homer/bin/findMotifsGenome.pl",
  result.dir = paste0("/media/RND/HDD-7/RZiffra/rziffra/allsnap/primary/AllPrimary_",levels(x.sp@cluster)[i],"_HOMER"),
  num.cores=16,
  genome = 'hg38',
  motif.length = 10,
  scan.size = 300,
  optimize.count = 25,
  local.background = FALSE,
  only.known = TRUE,
  only.denovo = FALSE,
  fdr.num = 1,
  cache = 100,
  overwrite = TRUE,
  keep.minimal = FALSE
);
}

library(chromVAR);
library(motifmatchr);
library(SummarizedExperiment);
library(BSgenome.Hsapiens.UCSC.hg38);
x.sp = makeBinary(x.sp, "pmat");
x.sp@mmat = runChromVAR(
  obj=x.sp,
  input.mat="pmat",
  genome=BSgenome.Hsapiens.UCSC.hg38,
  min.count=10,
  species="Homo sapiens"
);

motif_i = grep("GSX1",colnames(x.sp@mmat), value = T)[1];
dat = data.frame(x=x.sp@metaData[,"cluster"], y=x.sp@mmat[,motif_i]);
p1 <- ggplot(dat, aes(x=x, y=y, fill=x)) + 
  theme_classic() +
  geom_violin() + 
  xlab("cluster") +
  ylab("motif enrichment") +
  ylim(0,0.01) +
  ggtitle(motif_i) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.ticks.x=element_blank(),
    legend.position = "none"
  );
p1

par(mfrow=c(4,4))
for(i in 1:length(motif_idy)){
png(paste0(colnames(x.sp@mmat)[motif_idy[i]],".png"),res=300,height=3000,width=3000)
plotFeatureSingle(
  obj=x.sp,
  feature.value=x.sp@mmat[, motif_idy[8]],
  method="umap", 
  main=colnames(x.sp@mmat)[motif_idy[8]],
  point.size=1, 
  point.shape=19, 
  down.sample=30000,
  quantiles=c(0.05, 0.997))
dev.off()
}

library(rGREAT)
DARs = findDAR(
  obj=x.sp,
  input.mat="pmat",
  cluster.pos=1,
  cluster.neg=2,
  cluster.neg.method="knn",
  test.method="exactTest",
  bcv=0.2, #0.4 for human, 0.1 for mouse
  seed.use=10
);
DARs$FDR = p.adjust(DARs$PValue, method="BH");
idy = which(DARs$FDR < 5e-2 & DARs$logFC > 0);
job = submitGreatJob(
  gr                    = x.sp@peak[idy],
  bg                    = NULL,
  species               = "hg38",
  includeCuratedRegDoms = TRUE,
  rule                  = "basalPlusExt",
  adv_upstream          = 5.0,
  adv_downstream        = 1.0,
  adv_span              = 1000.0,
  adv_twoDistance       = 1000.0,
  adv_oneDistance       = 1000.0,
  request_interval = 300,
  max_tries = 10,
  version = "default",
  base_url = "http://great.stanford.edu/public/cgi-bin"
);
job

tb = getEnrichmentTables(job);
names(tb);
GBP = tb[["GO Biological Process"]];
head(GBP[order(GBP[,"Binom_Adjp_BH"]),],50);

pbmc.rna <- readRDS("/media/RND/HDD-7/RZiffra/rziffra/scRNA/GW17cortex/outs/filtered_feature_bc_matrix/GW17_Seurat.rds")
pbmc.rna$tech = "rna";
variable.genes = VariableFeatures(object = pbmc.rna)

pbmc <- snapToSeurat(
  obj=x.sp, 
  eigs.dims=1:32, 
  norm=TRUE,
  scale=TRUE
);

transfer.anchors <- FindTransferAnchors(
  reference = pbmc.rna, 
  query = pbmc.atac, 
  features = variable.genes[which(variable.genes %in% x.sp@gmat@Dimnames[[2]])], 
  reference.assay = "SCT", 
  query.assay = "ACTIVITY", 
  reduction = "cca"
);

celltype.predictions <- TransferData(
  anchorset = transfer.anchors, 
  refdata = pbmc.rna$seurat_clusters,
  weight.reduction = pbmc.atac[["SnapATAC"]],
  dims = 1:20
);

pbmc.atac <- AddMetaData(pbmc.atac, metadata = celltype.predictions)
pbmc.atac <- RunUMAP(pbmc.atac, reduction = "SnapATAC", dims = 1:20)
pbmc.atac.filtered <- subset(pbmc.atac, subset = prediction.score.max > 0.4)
pbmc.atac.filtered$predicted.id <- factor(pbmc.atac.filtered$predicted.id, levels = levels(pbmc.rna))

p1 <- DimPlot(pbmc.atac.filtered, group.by = "predicted.id", reduction = "umap", label = TRUE, repel = TRUE) + ggtitle("scATAC-seq cells") + 
  NoLegend() + scale_colour_hue(drop = FALSE)
p2 <- DimPlot(pbmc.rna, group.by = "seurat_clusters", label = TRUE, repel = TRUE) + ggtitle("scRNA-seq cells") + 
  NoLegend()
CombinePlots(plots = list(p1, p2))

refdata <- GetAssayData(
  object = pbmc.rna, 
  assay = "SCT", 
  slot = "data"
)
refdata <- refdata[variable.genes[which(variable.genes %in% x.sp@gmat@Dimnames[[2]])],]

imputation <- TransferData(
  anchorset = transfer.anchors, 
  refdata = refdata, 
  weight.reduction = pbmc.atac[["SnapATAC"]], 
  dims = 1:20
);

pbmc.atac[["SCT"]] <- imputation
DefaultAssay(pbmc.atac) <- "ATAC"
coembed <- merge(x = pbmc.rna, y = pbmc.atac)

coembed <- ScaleData(coembed, features = variable.genes[which(variable.genes %in% x.sp@gmat@Dimnames[[2]])], do.scale = FALSE)
coembed <- RunPCA(coembed, features = variable.genes[which(variable.genes %in% x.sp@gmat@Dimnames[[2]])], verbose = FALSE)
coembed <- RunUMAP(coembed, dims = 1:30)
coembed$seurat_clusters <- ifelse(!is.na(coembed$seurat_clusters), coembed$seurat_clusters, coembed$predicted.id)

p1 <- DimPlot(coembed, group.by = "tech")
p2 <- DimPlot(coembed, group.by = "seurat_clusters", label = TRUE, repel = TRUE)
CombinePlots(list(p1, p2))

x.sp@metaData$predicted.id = celltype.predictions$predicted.id;
x.sp@metaData$predict.max.score = apply(celltype.predictions[,-1], 1, max);
x.sp@cluster = as.factor(x.sp@metaData$predicted.id);

x.sp@gmat = t(imputation@data);
rm(imputation); # free memory
rm(refdata);    # free memory
rm(pbmc.rna);   # free memory
rm(pbmc.atac); # free memory

hist(
  x.sp@metaData$predict.max.score, 
  xlab="prediction score", 
  col="lightblue", 
  xlim=c(0, 1),
  main="PBMC 10X"
);
abline(v=0.5, col="red", lwd=2, lty=2);
table(x.sp@metaData$predict.max.score > 0.5);

x.sp = x.sp[x.sp@metaData$predict.max.score > 0.5,];
x.sp

plotViz(
  obj=x.sp,
  method="umap", 
  main="PBMC 10X",
  point.color=x.sp@metaData[,"predicted.id"], 
  point.size=0.5, 
  point.shape=19, 
  text.add=TRUE,
  text.size=1,
  text.color="black",
  down.sample=10000,
  legend.add=FALSE
)

x.sp = x.sp[which(x.sp@cluster %in% c(1,2,6,7,8,9,11,13,14,15,16,17,18,19,20,21))]
area_vector <-sapply(strsplit(x.sp@sample,"_"), `[`, 1)
runMACS(
  obj=x.sp[which(area_vector=="V1"),],
  output.prefix="PFCpeaks_EN",
  path.to.snaptools="/home/rziffra/.local/bin/snaptools",
  path.to.macs="/usr/local/bin/macs2",
  gsize="hs", # mm, hs, etc
  buffer.size=500,
  num.cores=8,
  macs.options="--nomodel --shift -37 --ext 73 --qval 5e-2 -B --SPMR",
  tmp.folder="/media/RND/HDD-7/RZiffra/rziffra/allsnap/primary/temp/"
)

