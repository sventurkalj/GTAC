# Projecting the query GTAC data to a reference scATAC-seq dataset; in our case, we projected AML data onto a healthy hematopoiesis reference

library(Seurat)
library(ggplot2)
library(patchwork)
library(Signac)
library(ArchR)
library(GenomicRanges)
library(GenomeInfoDb)
library(BSgenome.Hsapiens.UCSC.hg19)

set.seed(42)

# LOADING REFERENCE DATA

##### CREATING SEURAT OBJECT
ref   <- readRDS(file = "data/scATAC-Healthy-Hematopoiesis-191120.rds")

PEAKS_START <- c(ref@rowRanges@ranges@start)
PEAKS_END   <- c(ref@rowRanges@ranges@start)+500
PEAKS_CHR   <- c(as.vector((ref@rowRanges@seqnames)))

counts <- ref@assays$data$counts
colnames(counts) <- c(rownames(ref@colData))
rownames(counts) <- paste(PEAKS_CHR, paste(PEAKS_START, PEAKS_END, sep = "-"), sep = ":")

remove(PEAKS_START, PEAKS_END, PEAKS_CHR)

ranger <- StringToGRanges(regions = rownames(counts), sep = c(":", "-"))

chrom_assay <- CreateChromatinAssay(
  data = counts,
  min.cells = 0,
  min.features = 0,
  max.cells = NULL,
  ranges = ranger,
  motifs = NULL,
  fragments = NULL,
  genome = 'hg19',
  annotation = NULL,
  bias = NULL,
  positionEnrichment = NULL,
  sep = c(":", "-"),
  validate.fragments = NULL,
  verbose = TRUE
)

remove(counts, ranger)

new.ref <- CreateSeuratObject(counts=chrom_assay, project = "peaks", assay = "ATAC", row.names = colnames)

remove(chrom_assay)

df <- data.frame(rownames(ref@colData), ref@colData$BioClassification)
df$rownames.ref.colData <- gsub('_','-',df$rownames.ref.colData)
rownames(df) <- rownames(new.ref@meta.data)

new.ref <- AddMetaData(object = new.ref, metadata = df)

new.ref$input <- 'ref'

remove(df, ref)



##### LOADING MY DATA

query <- readRDS(file = "data/AML_only_hg19/Save-ArchR-Project.rds")

PEAKS_START <- c(query@peakSet@ranges@start)
PEAKS_END   <- c(query@peakSet@ranges@start)+500
PEAKS_CHR   <- c(as.vector((query@peakSet@seqnames)))

archr <- loadArchRProject(path = "data/AML_only_hg19", force = TRUE, showLogo = TRUE)

archr <- getMatrixFromProject(archr, useMatrix='PeakMatrix')
counts <- archr@assays
counts <- counts@data@listData[["PeakMatrix"]]

remove(PEAKS_START, PEAKS_END, PEAKS_CHR)

ranger <- archr@rowRanges

chrom_assay <- CreateChromatinAssay(
  data = counts,
  min.cells = 0,
  min.features = 0,
  max.cells = NULL,
  ranges = ranger,
  motifs = NULL,
  fragments = NULL,
  genome = 'hg19',
  annotation = NULL,
  bias = NULL,
  positionEnrichment = NULL,
  sep = c(":", "-"),
  validate.fragments = NULL,
  verbose = TRUE
)

remove(counts, ranger, archr)

new.query <- CreateSeuratObject(counts=chrom_assay, project = "peaks", assay = "ATAC", row.names = colnames)

remove(chrom_assay)

df <- data.frame(query$cellNames, query$clone)
df$query.cellNames <- gsub('_','-',df$query.cellNames)
rownames(df) <- rownames(new.query@meta.data)

new.query <- AddMetaData(object = new.query, metadata = df)

new.query$input <- 'query'

remove(query, df)



##### INTERSECTING PEAKS AND ADD QUERY UMAP

intersecting.regions <- GetIntersectingFeatures(
  object.1 = new.ref,
  object.2 = new.query,
  assay.1 = "ATAC",
  assay.2 = "ATAC",
  distance = 0
)

intersecting.regions.ref <- rownames(new.ref)[c(intersecting.regions[[1]])]
intersecting.regions.query <- rownames(new.query)[c(intersecting.regions[[2]])]

remove(intersecting.regions)

cells <- gsub("AML#", "", colnames(new.query))

new.query <- RenameCells(
    new.query[['ATAC']],
    new.names = cells
)

fragments <- CreateFragmentObject("data/AML_hg19_correct_sorted.tsv.gz", cells = cells)

remove(cells)

fragments <- FeatureMatrix(
  fragments = fragments,
  features = StringToGRanges(intersecting.regions.ref),
  cells = colnames(new.query)
)

new_chrom_assay <- CreateChromatinAssay(counts = fragments, genome = 'hg19')

remove(fragments)

new.query <- CreateSeuratObject(counts=new_chrom_assay, project = "peaks", assay = "ATAC")

remove(new_chrom_assay)

peaks.use <- intersecting.regions.ref #find subset of it

remove(intersecting.regions.ref, intersecting.regions.query)

query <- readRDS(file = "data/AML_only_hg19/Save-ArchR-Project.rds")

df <- data.frame(query$cellNames, query$clone)
df$query.cellNames <- gsub('_','-',df$query.cellNames)
rownames(df) <- rownames(new.query@meta.data)

new.query <- AddMetaData(object = new.query, metadata = df)

UMAP_1 <- query@embeddings$UMAP_LSI_peaks$df$`IterativeLSI_peaks#UMAP_Dimension_1`
UMAP_2 <- query@embeddings$UMAP_LSI_peaks$df$`IterativeLSI_peaks#UMAP_Dimension_2`
umap <- cbind(UMAP_1, UMAP_2)
new.query[["query.umap"]] <- CreateDimReducObject(embeddings = umap, key = "UMAP_", assay = DefaultAssay(new.query))
rownames(new.query@reductions$query.umap@cell.embeddings) <- colnames(new.query)

new.query$input <- 'query'

remove(umap, UMAP_1, UMAP_2, df, query)

new.ref <- new.ref[rownames(new.ref) %in% peaks.use, ]
new.query <- new.query[rownames(new.query) %in% peaks.use, ]

peaks.use <- rownames(new.query[new.query[['ATAC']]@meta.features$count>0])

new.ref <- new.ref[rownames(new.ref) %in% rownames(new.query[new.query[['ATAC']]@meta.features$count>0]), ]
new.query <- new.query[rownames(new.query) %in% rownames(new.query[new.query[['ATAC']]@meta.features$count>0]), ]



##### CALCULATING LSIs

d <- 30

new.query <- RunTFIDF(new.query)
new.query <- FindTopFeatures(new.query, min.cutoff = 311) #10% 311.4
new.query <- RunSVD(new.query, n = d, reduction.name = 'lsi', reduction.key = 'LSI_')

new.ref@assays$ATAC@counts <- new.ref[['ATAC']][]
new.ref@assays$ATAC@data <- RunTFIDF(new.ref[['ATAC']][])
new.ref <- FindTopFeatures(new.ref, min.cutoff = 3504) #10% 3503.8
new.ref <- RunSVD(new.ref, n = d, reduction.name = 'lsi', reduction.key = 'LSI_')

DefaultAssay(new.ref) <- "ATAC"
DefaultAssay(new.query) <- "ATAC"



##### SAVING rds FILES 1


saveRDS(new.ref, file = "rds/new.ref.1.rds")
saveRDS(new.query, file = "rds/new.query.1.rds")



##### CHECKPOINT (if files have been saved, you can start from here, removing #)

#library(Seurat)
#library(ggplot2)
#library(patchwork)
#library(Signac)
#library(ArchR)
#library(GenomicRanges)
#library(GenomeInfoDb)
#library(BSgenome.Hsapiens.UCSC.hg19)

#set.seed(42)

#new.ref <- readRDS(file = "rds/new.ref.1.rds")
#new.query <- readRDS(file = "rds/new.query.1.rds")

#d <- 30
#peaks.use <- rownames(new.ref)



##### ADDING UMAP MODEL

umapManifold <- uwot::load_uwot("data/scATAC-Hematopoiesis-UMAP-model.190505.uwot.tar")
rownames(umapManifold$embedding) <- colnames(new.ref)
new.ref[["umap"]] <- CreateDimReducObject(embeddings = umapManifold$embedding, key = "UMAP_", assay = DefaultAssay(new.ref), misc=list(umapManifold))
names(new.ref@reductions$umap@misc) <- 'model'

remove(umapManifold)



##### REMOVING UNK CLUSTERS

Idents(object = new.ref) <- 'ref.colData.BioClassification'
new.ref <- subset(x = new.ref, idents = c("26_Unk", "14_Unk", "13_Unk"), invert = TRUE)

new.ref@reductions$umap@misc$model$embedding <- new.ref@reductions$umap@cell.embeddings



##### FINIDING ANCHORS

transfer.anchors <- FindTransferAnchors(
  reference = new.ref,
  query = new.query,
  reference.assay = 'ATAC',
  query.assay = 'ATAC',
  reduction = 'lsiproject',
  reference.reduction = 'lsi',
  dims = 2:d,
  features = peaks.use
)



##### MAPPING QUERY

new.query <- MapQuery(
  anchorset = transfer.anchors,
  query = new.query,
  reference = new.ref,
  refdata = list(predicted.celltype = "ref.colData.BioClassification"),
  reference.reduction = "lsi",
  reduction.model = "umap"
)



##### SAVING rds FILES 2

saveRDS(new.ref, file = "rds/new.ref.2.rds")
saveRDS(new.query, file = "rds/new.query.2.rds")



##### CONVERTING FROM Seurat TO AnnData via h5Seurat

library(SeuratDisk)

SaveH5Seurat(new.query, filename = "h5ad/query.h5Seurat", overwrite = TRUE)
Convert("h5ad/query.h5Seurat", dest = "h5ad")

new.ref.2 <- new.ref
new.ref.2[["umap"]] <- CreateDimReducObject(embeddings = new.ref@reductions$umap@cell.embeddings, key = "UMAP_", assay = DefaultAssay(new.ref))

SaveH5Seurat(new.ref.2, filename = "h5ad/ref.h5Seurat", overwrite = TRUE)
Convert("h5ad/ref.h5Seurat", dest = "h5ad")


##### SAVING INTERSECTING PEAKS

write.table(peaks.use, "intersecting_peaks.csv", sep="\t")



