############################################
# ATAC batch correction
# After running LSI/LDA/etc. in SnapATAC
mtx <- x.sp@smat@dmat
rownames(mtx) <- x.sp@barcode
colnames(x.sp) <- 1:50 # B/c seurat needs names; use however many dims needed
seuratObj <- CreateSeuratObject(t(mtx))
seuratObj@meta.data$sample <- x.sp@metaData$sample # Batch info

seurat.list <- SplitObject(seuratObj, split.by='sample')
sce.list <- lapply(seurat.list, function(i){
  SingleCellExperiment(assays=list(LSI = i@assays$RNA@counts))
})
sca <- scAlignCreateObject(sce.list)
rm(sce.list);gc()

# Generally speaking 15k is enough to converge, but if there are a lot of batch effect 20k+ might be required
# Architecture size doesn't really matter, small is good enough most of the time. Increase to 'large' if it's not correcting enough
# If it's overaligned, turn norm=F and turn batch.norm.layer=T. This will reduce correction
# num.dim should generally be lower than your initial dim (which is 50 here) unless you are using CCA.
sca = scAlignMulti(sca,
                   options=scAlignOptions(steps=15000,norm=TRUE, batch.norm.layer=F, early.stop=FALSE, architecture='small', num.dim=32),
                   encoder.data="LSI",
                   supervised='none',
                   run.encoder=TRUE,
                   run.decoder=F,
                   device="CPU") # Or GPU, but CPU is faster, this will use all cores available
batch.corrected <- reducedDim(sca, "ALIGNED-GENE") # ALIGNED-GENE is currently hardcoded into scAlign
rownames(batch.corrected) <- colnames(seuratObj) 
seuratObj[['sca']] <- CreateDimReducObject(batch.corrected, key="SCA_")
seuratObj <- RunUMAP(seuratObj, dims=1:32, reduction='sca') # Seurat UMAP parameters are ideal for most single cell datasets, SnapATAC is no bueno
umap <- Embeddings(seuratObj, reduction='sca') # Pull UMAP coords

barcodes <- make.unique(x.sp@barcode) # SnapATAC for some reason doesn't actually make barcodes unique like Seurat, skip if it's already made unique
x.sp@umap <- umap[barcodes,] # This will order the umap coord with the proper cells
x.sp <- RunLSI(x.sp, logTF=T, pc.num=32) # Rerun num.dim used in scAlign since SnapATAC won't allow different dims
x.sp@smat@dmat <- batch.corrected[barcodes,]

##########################################
# RNA-ATAC integration
# Run SnapATAC all the way to generating gene activity scores and specifically use the raw RPM values
# Do not use MAGIC
atac <- CreateSeuratObject(gmat)
atac <- NormalizeData(atac)
# Process RNA using SCTransform if UMI based or the usual NormalizeData
rna <- SCTransform(rna)
feats <- SelectIntegrationFeatures(list(rna,atac), nfeatures = 3000) 
rna <- GetResidual(rna, feats)
VariableFeatures(atac) <- feats
atac <- ScaleData(atac)
sce.list <- lapply(list(rna,atac), function(i){
  if(DefaultAssay(i) == "RNA"){
    SingleCellExperiment(assays=list(scale.data = i@assays$RNA@scale.data[feats,]))
  }
  else{
    SingleCellExperiment(assays=list(scale.data = i@assays$SCT@scale.data[feats,]))
  }
})
sca <- scAlignCreateObject(sce.list, cca.reduce=T,ccs.compute=10) # Num of CCs chosen matters a lot, I usually use 10-15 depending on how similar the datasets are
sca = scAlignMulti(sca,
                   options=scAlignOptions(steps=15000,
                                          log.every=5000,
                                          batch.size=300, # 300 is for larger datasets, default is 100. Speeds up training with higher values
                                          norm=TRUE,
                                          batch.norm.layer=FALSE,
                                          architecture="small", 
                                          num.dim=20), # Use small-ish values.    
                   encoder.data="MultiCCA",
                   supervised='none',
                   run.encoder=TRUE,
                   device="CPU")
batch.corrected <- reducedDim(sca, "ALIGNED-MultiCCA") 
merged <- merge(rna, atac) # maintain order of the objects used to create scAlign obj
rownames(batch.corrected) <- colnames(merged)
merged[['sca']] <- CreateDimReducObject(batch.corrected, key="SCA_")
merged <- RunUMAP(merged, reduction='sca')
merged <- FindNeighbors(merged, reduction='sca', dims=1:32)
